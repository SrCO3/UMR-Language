import json
from ..chemistry import FormalError,Unsupported,atoms_by_id,bond_orders,electron_count
from ..v05.ledger import validate_ledger as ledger05,compare_views,signature
from ..v04.ledger import validate_view as validate04
from .elements import MAIN,METAL_GROUP,coordination,local_profile,oxidation,electronegativities


def symbols(view):return {a.atom_id for a in view.atoms if a.element in {'*','R'}}


def project_symbols(view):
    """Private boundary caps for local checking; never exported as chemistry."""
    ids=symbols(view);shadow=view.model_copy(deep=True)
    if ids and view.annotations.get('symbolic_mode')!='unchanged_boundary':raise Unsupported('Symbolic execution must be explicitly enabled')
    for atom in shadow.atoms:
        if atom.atom_id not in ids:continue
        edges=[b for b in view.bonds if atom.atom_id in {b.a,b.b}]
        if len(edges)!=1 or edges[0].order!=1 or atom.hydrogens or atom.isotope:
            raise Unsupported('Opaque R/n scope: one unchanged single-bond boundary, no known internal H/isotope')
        if atom.formal_charge is not None or atom.nonbonding_electrons is not None or atom.radical_electrons is not None:
            raise FormalError('Symbolic fragment occupancy must remain unknown')
        atom.element='F';atom.formal_charge=0;atom.nonbonding_electrons=6;atom.radical_electrons=0
    return shadow


def restore_symbols(original,target):
    ids=symbols(original);old=atoms_by_id(original)
    a={k:v for k,v in bond_orders(original).items() if set(k)&ids};b={k:v for k,v in bond_orders(target).items() if set(k)&ids}
    if a!=b:raise FormalError('Opaque fragment boundary changed')
    target.atoms=[old[x.atom_id].model_copy(deep=True) if x.atom_id in ids else x for x in target.atoms]
    return target


def validate_ledger(view):
    shadow=project_symbols(view)
    # Radical Sn/Se occupy their own finite profiles, not the old organic checker.
    report=ledger05(shadow,_profile=local_profile,_coordination=coordination,
                    _main=MAIN,_metals=METAL_GROUP,_oxidation=lambda v,i:oxidation(v,i)) if not any(a.element in {'Sn','Se'} and a.radical_electrons for a in shadow.atoms) else _radical_extended_ledger(shadow)
    report['rule_version']='umr-elements/0.7'
    report['validation_scope']='symbolic_unchanged_boundary' if symbols(view) else 'complete_material_ledger'
    if symbols(view):
        report['opaque_symbols']=sorted(symbols(view));report['full_electron_total_verified']=False
        report['known_electrons_plus_opaque_terms']=report['electrons']-6*len(symbols(view))
        report['electrons']=None
    table=electronegativities()
    report['oxidation_analysis']=[{'atom_id':a.atom_id,'oxidation_state':oxidation(view,a.atom_id),
                                 'scale':'Allen','status':'derived' if oxidation(view,a.atom_id) is not None else 'unknown',
                                 'parameter_record_ids':sorted({table[x.element]['record_id'] for x in view.atoms if x.element in table}|({table['H']['record_id']} if 'H' in table and any(x.hydrogens for x in view.atoms) else set()))}
                                for a in view.atoms if a.element not in METAL_GROUP]
    return report


def _radical_extended_ledger(view):
    checks=[];shadow=view.model_copy(deep=True)
    for a in shadow.atoms:
        if a.element in {'Sn','Se'} and a.radical_electrons:
            checks.append(local_profile(view,atoms_by_id(view)[a.atom_id]))
            # Element-local equivalent organic valence solely for legacy validation.
            a.element='C' if a.element=='Sn' else 'O'
    result=ledger05(shadow,_profile=local_profile,_coordination=coordination,_main=MAIN,
                    _metals=METAL_GROUP,_oxidation=oxidation)
    result['atom_checks']=[x for x in result['atom_checks'] if x['atom_id'] not in {c['atom_id'] for c in checks}]+checks
    return result


def validate_view(view,state=None):
    report=validate_ledger(view);shadow=project_symbols(view)
    unknown={a.atom_id for a in shadow.atoms if a.radical_electrons is None}
    if unknown:
        if any(atoms_by_id(shadow)[i].element not in METAL_GROUP for i in unknown):raise Unsupported('Unknown active main-group spin')
        if view.electronic_state.character!='unknown' or view.electronic_state.coupling!='unknown' or (state and state.multiplicity is not None):
            raise Unsupported('Unknown metal spin cannot establish a definite global electronic state')
        for a in shadow.atoms:
            if a.atom_id in unknown:a.radical_electrons=0
    checking={**report,'total_charge':sum(a.formal_charge for a in shadow.atoms),
              'unpaired_electrons':sum(a.radical_electrons for a in shadow.atoms)}
    validate04(shadow,state,_ledger=lambda v:checking)
    report['spin_status']='unknown_not_inferred' if unknown else 'declared_occupation_only'
    return report
