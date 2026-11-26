"""Local stereochemical reassessment; uncertainty is never an inversion proof."""
RULE_VERSION="umr-stereo-review/0.1"


def stereo_assertion(value):
    # Reasons/evidence on an unknown assignment are annotations, not a geometry.
    if value.get("status")!="specified":return None
    return {k:v for k,v in value.items() if k not in {"reason","evidence","note"}}


def reassess_centers(source,target,touched,*,resonance=False):
    original={a.atom_id:a for a in source.atoms};notes=[]
    for atom in target.atoms:
        if atom.atom_id not in touched:continue
        before=original[atom.atom_id]
        labile=atom.element=="C" and (atom.radical_electrons or atom.formal_charge in {-1,1} or before.radical_electrons or before.formal_charge in {-1,1})
        changed=any(b.order>1 and atom.atom_id in {b.a,b.b} for b in target.bonds)
        if resonance and not (labile or changed):continue
        # A deliberately supplied stability/geometry constraint must not be erased.
        stable=bool(before.stereo.get("configurational_stability_evidence"))
        if resonance and stable:continue
        reason="该原子参与电子转移，需要重新判断立体构型；缺少三维/立体证据"
        if labile:reason+="；自由基或带电中心的构型可能不确定，不能仅凭二维图赋值"
        if resonance:reason+="；此处是同一电子态的共振表示，不代表真实翻转"
        atom.stereo={"status":"unknown","reason":reason}
        notes.append({"atom_id":atom.atom_id,"status":"undetermined","before":before.stereo,
            "reason":reason,"rule_version":RULE_VERSION,"explicit_stable_center_not_generalized":True})
    if notes:target.annotations["stereo_reassessment"]=notes
    return notes


def check_resonance_stereo(source,target):
    from .chemistry import FormalError
    before={a.atom_id:a for a in source.atoms}
    after={a.atom_id:a for a in target.atoms}
    for i,a in before.items():
        old,new=stereo_assertion(a.stereo),stereo_assertion(after[i].stereo)
        if old==new:continue
        same_neighbors=lambda v:sorted((min(b.a,b.b),max(b.a,b.b),b.order) for b in v.bonds if i in {b.a,b.b})
        local=(same_neighbors(source)!=same_neighbors(target) or
            (a.radical_electrons,a.formal_charge,a.nonbonding_electrons)!=(after[i].radical_electrons,after[i].formal_charge,after[i].nonbonding_electrons))
        # Withdrawing an unsupported assignment at a local labile center is
        # uncertainty, not inversion. A new specified geometry is not authorized.
        if old is not None and new is None and local and not a.stereo.get("configurational_stability_evidence"):continue
        raise FormalError("Resonance cannot invert atom stereochemistry at an established center")
