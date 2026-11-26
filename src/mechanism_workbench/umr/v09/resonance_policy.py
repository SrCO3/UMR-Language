"""Contributor selection, independent of formal validity and real ionization."""
from dataclasses import dataclass
from math import isfinite

VERSION='umr-resonance-selection/0.1'

@dataclass(frozen=True)
class ResonancePolicy:
    suppress_neutral_cc_separation: bool=True
    min_site_electronegativity_difference: float=0.25
    def __post_init__(self):
        if not isfinite(self.min_site_electronegativity_difference) or self.min_site_electronegativity_difference<=0:
            raise ValueError('Positive finite site-electronegativity threshold required')

def decision(source,target,policy):
    if not policy.suppress_neutral_cc_separation:return None
    aa={a.atom_id:a for a in source.atoms};bb={a.atom_id:a for a in target.atoms}
    positive=[i for i,a in aa.items() if a.element=='C' and a.formal_charge==0 and bb[i].formal_charge>0]
    negative=[i for i,a in aa.items() if a.element=='C' and a.formal_charge==0 and bb[i].formal_charge<0]
    if not positive or not negative:return None
    adjacency={i:set() for i in aa}
    for b in source.bonds:adjacency[b.a].add(b.b);adjacency[b.b].add(b.a)
    component=set(positive);todo=list(positive)
    while todo:
        for j in adjacency[todo.pop()]-component:component.add(j);todo.append(j)
    if any(aa[i].formal_charge for i in component):return None
    raw=source.annotations.get('site_electronegativity',{})
    report={'policy_version':VERSION,'positive_carbon_ids':positive,'negative_carbon_ids':negative,
        'threshold':policy.min_site_electronegativity_difference,
        'threshold_status':'engineering_selection_default_not_a_universal_chemical_constant',
        'source_view_id':source.view_id,'formal_validity':'independent'}
    if not isinstance(raw,dict) or not all(raw.get(k) for k in ('scale','model_version','source','values')):
        return {**report,'selected':False,'reason':'equal_elemental_C_electronegativity_and_no_supported_site_difference'}
    if raw['scale']!='Pauling-equivalent' or raw.get('calibrated') is not True:
        return {**report,'selected':False,'reason':'site_scale_or_calibration_not_supported'}
    values=raw['values'];diffs=[]
    if not isinstance(values,dict):
        return {**report,'selected':False,'reason':'site_values_must_be_atom_keyed_mapping'}
    try:
        for p in positive:
            for n in negative:
                if n not in component:continue
                if isinstance(values[str(n)],bool) or isinstance(values[str(p)],bool):raise ValueError()
                x,y=float(values[str(n)]),float(values[str(p)])
                if not isfinite(x) or not isfinite(y):raise ValueError()
                diffs.append(x-y)
    except (KeyError,ValueError,TypeError):
        return {**report,'selected':False,'reason':'site_values_missing_or_nonfinite'}
    return {**report,'selected':bool(diffs) and min(diffs)>policy.min_site_electronegativity_difference,
        'signed_differences':diffs,'source':raw['source'],'model_version':raw['model_version'],
        'reason':'negative_carbon_must_have_sufficiently_higher_supported_site_electronegativity'}
