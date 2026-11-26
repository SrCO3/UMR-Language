import json
from ..ast import Bond
from ..chemistry import FormalError, Unsupported, atoms_by_id, bond_valences, electron_count
from ..v04.ast import Atom, LewisView as View04
from ..v04.ledger import validate_ledger as ledger04, validate_view as validate04, signature as signature04, compare_views as compare04
from .elements import MAIN, METAL_GROUP, charge, local_profile, oxo_oxidation
from .coordination import checks as coordination_checks


def extension_signature(view):
    return json.dumps({k:getattr(view,k) if k=='element_rule_version' else [x.model_dump(mode='json') for x in getattr(view,k)]
        for k in ('element_rule_version','element_declarations','lewis_adducts')},sort_keys=True)


def signature(view):return signature04(view),extension_signature(view)


def compare_views(actual,expected):
    compare04(actual,expected)
    if extension_signature(actual)!=extension_signature(expected):raise FormalError('Element/coordination declarations differ from asserted target')


def validate_ledger(view,*,_profile=local_profile,_coordination=coordination_checks,
                    _main=MAIN,_metals=METAL_GROUP,_oxidation=oxo_oxidation,_local_atoms=()):
    MAIN, METAL_GROUP = _main, _metals
    ids={a.atom_id for a in view.atoms};edges=[tuple(sorted((b.a,b.b))) for b in view.bonds]
    if len(ids)!=len(view.atoms) or len(set(edges))!=len(edges) or any(b.a==b.b or not {b.a,b.b}<=ids or b.order not in {1,2,3} for b in view.bonds):
        raise FormalError('Atom IDs/bond references must be unique, complete and typed')
    if any(d.atom_id not in ids for d in view.element_declarations) or any(not {a.donor_atom_id,a.acceptor_atom_id}<=ids for a in view.lewis_adducts):
        raise FormalError('Element/adduct reference is absent')
    if len({m.atom_id for m in view.metal_centers})!=len(view.metal_centers):raise FormalError('Duplicate metal declaration')
    vals=bond_valences(view);special=set();reports=[]
    for a in view.atoms:
        if a.atom_id in _local_atoms or a.element in METAL_GROUP or (a.element in MAIN and not a.radical_electrons):
            reports.append(_profile(view,a));special.add(a.atom_id)
    reserved,counts=_coordination(view)
    # Reuse the complete v0.4 organic checks on a PRIVATE capped fragment. These
    # cap atoms are never exported or used in the whole-system electron count.
    ordinary=[a.model_copy(deep=True) for a in view.atoms if a.atom_id not in special]
    bonds=[b.model_copy(deep=True) for b in view.bonds if not {b.a,b.b}&special]
    next_id=max(a.atom_id for a in view.atoms)+1;caps=[]
    for a in ordinary:
        lost=sum(b.order for b in view.bonds if a.atom_id in {b.a,b.b} and {b.a,b.b}&special)
        if a.element=='H' and lost:
            if lost!=1:raise FormalError('Hydrogen may have only one single bond')
            ordinary_cap=Atom(atom_id=next_id,element='H',formal_charge=0,nonbonding_electrons=0)
            bonds.append(Bond(a=a.atom_id,b=next_id,order=1));next_id+=1
            # append after iterating the original list
            caps.append(ordinary_cap)
        else:a.hydrogens+=lost
    if ordinary:
        shadow=View04(state_id=view.state_id,view_id=view.view_id,atoms=ordinary+caps,bonds=bonds,
                      reactive_centers=view.reactive_centers,aryne_regions=view.aryne_regions,electronic_state=view.electronic_state)
        reports+=ledger04(shadow)['atom_checks']
    byid={r['atom_id']:r for r in reports};atoms=atoms_by_id(view)
    for d in view.element_declarations:
        a=atoms[d.atom_id]
        cn=sum(d.atom_id in {b.a,b.b} for b in view.bonds)+a.hydrogens
        cn+=sum(x.acceptor_atom_id==d.atom_id or x.donor_atom_id==d.atom_id for x in view.lewis_adducts)
        cn+=sum(d.atom_id in b.donor_atom_ids for b in view.metal_bindings if b.mode=='dative')
        cn+=sum(len(b.donor_atom_ids) for b in view.metal_bindings if b.metal_atom_id==d.atom_id and b.mode!='sigma_covalent')
        if d.coordination_number is not None and d.coordination_number!=cn:raise FormalError('Declared coordination number differs from graph/bindings')
        os=next((c['oxidation_state'] for c in counts if c['metal_atom_id']==d.atom_id),None) if a.element in METAL_GROUP else _oxidation(view,d.atom_id)
        if d.oxidation_state is not None:
            if os is None:raise Unsupported('Oxidation state cannot be uniquely checked in this environment')
            if d.oxidation_state!=os:raise FormalError('Oxidation state differs from ionic approximation, not from formal charge')
        if d.atom_id in byid:byid[d.atom_id].update(oxidation_state=os,declared_valence_category=d.valence_category,coordination_number=cn)
    return {'total_charge':sum(a.formal_charge for a in view.atoms),'electrons':electron_count(view),
            'unpaired_electrons':None if any(a.radical_electrons is None for a in view.atoms) else sum(a.radical_electrons for a in view.atoms),
            'atom_checks':reports,'metal_counts':counts,'reserved_donor_electrons':dict(reserved),
            'rule_version':'umr-elements/0.5','chemical_feasibility':'unknown'}


def validate_view(view,state=None,*,_profile=local_profile):
    report=validate_ledger(view,_profile=_profile)
    unknown={a.atom_id for a in view.atoms if a.radical_electrons is None}
    if unknown:
        if any(atoms_by_id(view)[i].element not in METAL_GROUP for i in unknown):raise Unsupported('Unknown nonmetal occupation')
        if view.electronic_state.character!='unknown' or view.electronic_state.coupling!='unknown' or (state and (state.multiplicity is not None or state.electronic_mode!='UNKNOWN')):
            raise Unsupported('Unresolved metal spin cannot validate a definite electronic-state declaration')
        shadow=view.model_copy(deep=True)
        for a in shadow.atoms:
            if a.atom_id in unknown:a.radical_electrons=0
        # Only run the legacy declaration guards on this private compatibility
        # view. Unknown occupations in the authoritative view remain unknown.
        validate04(shadow,state,_ledger=lambda v:{**report,'unpaired_electrons':sum(a.radical_electrons for a in shadow.atoms)})
        report['spin_status']='unknown_not_inferred'
    else:
        validate04(view,state,_ledger=lambda v:report)
    return report
