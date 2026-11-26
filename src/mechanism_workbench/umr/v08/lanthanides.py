"""Finite ionic La coordination and declared Ce/Sm outer-sphere SET accounting.

No f-shell wavefunctions, PCET, Ce(II) generic reduction or ligand redox inference.
"""
from ..v05.ast import Ligand, MetalBinding
from ..chemistry import FormalError, Unsupported, atoms_by_id, bond_valences, identity, electron_count
from .ledger import validate_view


def execute(source,event):
    if not event.evidence or not event.conditions or event.actions is not None:
        raise FormalError('Finite lanthanide recipe requires evidence/conditions and compiles its own electron accounting')
    target=source.model_copy(deep=True);target.state_id=event.product_state_id;target.view_id=event.product_view_id
    before=atoms_by_id(source);after=atoms_by_id(target);name=event.element_template
    if name.startswith('LA_L_'):
        if event.event_type!='COORDINATION' or set(event.bindings)!={'metal','donor'}:raise FormalError('La roles/event type mismatch')
        m,d=event.bindings['metal'],event.bindings['donor']
        if m not in before or d not in before or (before[m].element,before[m].formal_charge)!=('La',3) or before[d].element not in {'O','N'}:
            raise Unsupported('La(III), neutral connected N/O donor profile only')
        existing=[b for b in source.metal_bindings if b.metal_atom_id==m and b.donor_atom_ids==[d]]
        if name=='LA_L_DISSOCIATION':
            if len(existing)!=1:raise FormalError('Missing unique La donor binding')
            binding=existing[0];target.metal_bindings=[b for b in target.metal_bindings if b.binding_id!=binding.binding_id]
            target.ligands=[l for l in target.ligands if l.ligand_id!=binding.ligand_id]
        else:
            if existing:raise FormalError('Repeated La donor coordination')
            members={d};stack=[d]
            while stack:
                i=stack.pop()
                for b in source.bonds:
                    if i not in {b.a,b.b}:continue
                    j=b.b if b.a==i else b.a
                    if j not in members:members.add(j);stack.append(j)
            if m in members or sum(before[i].formal_charge for i in members)!=0:raise Unsupported('La donor must be a neutral separate ligand')
            lid=f'la.ligand.{d}';bid=f'la.binding.{m}.{d}'
            target.ligands.append(Ligand(ligand_id=lid,atom_ids=sorted(members),donor_atom_ids=[d],ionic_class='L',ionic_charge=0))
            target.metal_bindings.append(MetalBinding(binding_id=bid,metal_atom_id=m,ligand_id=lid,donor_atom_ids=[d],mode='dative'))
        transfer={'kind':'donor_owned_coordination','metal':m,'donor':d,'global_electrons_added':0}
    else:
        if event.event_type!='SET' or set(event.bindings)!={'donor','acceptor'}:raise FormalError('SET needs donor/acceptor roles')
        d,a=event.bindings['donor'],event.bindings['acceptor']
        if d==a or not {d,a}<=before.keys():raise FormalError('Missing/distinct SET centers required')
        metal=a if name=='CEIV_ONE_ELECTRON_OXIDATION' else d
        organic=d if metal==a else a
        el,q=('Ce',4) if metal==a else ('Sm',2)
        if (before[metal].element,before[metal].formal_charge)!=(el,q) or bond_valences(source)[metal] or any(b.metal_atom_id==metal for b in source.metal_bindings):
            raise Unsupported('First SET is declared ionic outer-sphere CeIV/SmII; coordinated or noninnocent ligands need a separate recipe')
        atom=before[organic]
        if atom.element!='C' or bond_valences(source)[organic]!=3 or (atom.formal_charge,atom.nonbonding_electrons,atom.radical_electrons)!=((-1,2,0) if organic==d else (1,0,0)):
            raise Unsupported('Formal SET fixture scope is localized carbanion oxidation or carbocation reduction')
        if not event.conditions.get('electron_transfer_model')=='outer_sphere_formal':
            raise FormalError('Explicit outer-sphere model required; SmI2/water is not automatically isolated SET')
        after[d].nonbonding_electrons-=1;after[d].formal_charge+=1
        after[a].nonbonding_electrons+=1;after[a].formal_charge-=1
        after[organic].radical_electrons=1;after[metal].radical_electrons=None
        for m in target.metal_centers:
            if m.atom_id==metal:m.oxidation_state=after[metal].formal_charge;m.d_electrons=None;m.valence_electrons=None;m.spin_state='unknown'
        for dec in target.element_declarations:
            if dec.atom_id==metal:dec.oxidation_state=after[metal].formal_charge
        transfer={'kind':'single_electron','donor':d,'acceptor':a,'electrons':1,'spin_coupling':'unknown','pcet_executed':False}
    validate_view(target)
    if identity(source)!=identity(target) or electron_count(source)!=electron_count(target) or sum(a.formal_charge for a in source.atoms)!=sum(a.formal_charge for a in target.atoms):
        raise FormalError('Lanthanide recipe violates material/electron/charge balance')
    return target,{'event_id':event.event_id,'event_type':event.event_type,'execution_mode':'electron_actions',
        'typed_electron_transfer':transfer,'compiled_actions':[],'atomic_commit':True,'formal_status':'passed',
        'scope':'finite explicit ionic recipe; coordination/spin dynamics not computed',
        'chemical_feasibility':'unknown','human_approved':False}
