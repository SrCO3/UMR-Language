"""One-turnover receipts are derived from executed states, not model claims."""
import hashlib,json
from ..chemistry import FormalError,identity,electron_count,bond_orders
from .ledger import validate_view

def mapped_fragment(view,ids,mapping=None):
    ids=set(ids);mapping=mapping or {i:i for i in ids}
    stereo=lambda s:json.dumps(s if s.get('status')=='specified' else {'status':'unknown'},sort_keys=True)
    return {'atoms':sorted((mapping[a.atom_id],a.element,a.isotope,a.formal_charge,a.nonbonding_electrons,a.radical_electrons,a.hydrogens,stereo(a.stereo)) for a in view.atoms if a.atom_id in ids),
       'bonds':sorted((min(mapping[b.a],mapping[b.b]),max(mapping[b.a],mapping[b.b]),b.order,stereo(b.stereo)) for b in view.bonds if {b.a,b.b}<=ids),
       'coordination':sorted((mapping[b.metal_atom_id],tuple(sorted(mapping[i] for i in b.donor_atom_ids)),b.mode,b.hapticity) for b in view.metal_bindings if {b.metal_atom_id,*b.donor_atom_ids}<=ids)}
def validate_cycle(document,reports,cycle):
    events={e.event_id:e for e in document.events};views={v.view_id:v for v in document.views};run={r['event_id']:r for r in reports}
    allids=cycle.event_ids+cycle.activation_event_ids+cycle.termination_event_ids
    if len(set(cycle.event_ids))!=len(cycle.event_ids) or not set(allids)<=events.keys():raise FormalError('Missing/repeated cycle event reference')
    if not set(cycle.event_ids)<=run.keys():return {'status':'incomplete_cycle','compressible':False,'reason':'unexecuted turnover events','reward':0}
    if any(run[i].get('execution_mode','electron_actions')!='electron_actions' or run[i].get('depends_on_exceptions') or run[i].get('depends_on_net_events') or run[i].get('occupation_assumption') or run[i].get('mechanistic_resolution')=='grouped_redox_proton_ledger_not_asserted_elementary' for i in cycle.event_ids):
        return {'status':'not_strict_cycle','compressible':False,'reward':0}
    current=cycle.initial_view_id
    for i in cycle.event_ids:
        e=events[i]
        if e.reactant_view_id!=current:raise FormalError('Cycle events are not a contiguous executed path')
        current=e.product_view_id
    if current!=cycle.final_view_id:raise FormalError('Cycle final view is not actual executed endpoint')
    s,t=views[cycle.initial_view_id],views[cycle.final_view_id];validate_view(s);validate_view(t)
    from ..v07.ledger import symbols
    if symbols(s) or symbols(t):return {'status':'not_strict_cycle','compressible':False,'reason':'opaque inventory prevents complete cycle ledger','reward':0}
    if any(b.direction!='inside' for v in (s,t) for b in v.material_boundaries):return {'status':'incomplete_cycle','compressible':False,'reason':'external exchange closure not verified','reward':0}
    if s.electronic_state!=t.electronic_state:raise FormalError('Cycle electronic-state declarations differ')
    if identity(s)!=identity(t) or electron_count(s)!=electron_count(t) or sum(a.formal_charge for a in s.atoms)!=sum(a.formal_charge for a in t.atoms):raise FormalError('Cycle material/electron/charge mismatch')
    mapping=cycle.carrier_mapping
    if not mapping or len(set(mapping.values()))!=len(mapping):raise FormalError('Carrier role mapping must be nonempty bijection')
    if not set(mapping)<={a.atom_id for a in s.atoms} or not set(mapping.values())<={a.atom_id for a in t.atoms}:raise FormalError('Missing carrier atom')
    if any(mapping.get(i)!=i for i in cycle.persistent_core_ids):raise FormalError('Persistent catalytic core identity changed')
    if {i for i,j in mapping.items() if i!=j}-set(cycle.exchangeable_atom_ids):raise FormalError('Carrier exchange lacks explicit scope')
    reverse={v:k for k,v in mapping.items()}
    if mapped_fragment(s,mapping)!=mapped_fragment(t,mapping.values(),reverse):raise FormalError('Carrier not regenerated under declared role correspondence')
    groups=cycle.product_atom_groups;seen=set();changed=False
    for group in groups:
        ids=set(group)
        if not ids or len(ids)!=len(group) or seen&ids or not ids<={a.atom_id for a in t.atoms} or ids&set(mapping.values()):raise FormalError('Invalid/duplicate product inventory')
        seen|=ids
        if any(bool(b.a in ids)!=bool(b.b in ids) for b in t.bonds) or any(bool(b.metal_atom_id in ids)!=bool(set(b.donor_atom_ids)&ids) for b in t.metal_bindings):raise FormalError('Product has not been released from inventory/carrier')
        changed|=mapped_fragment(s,group)!=mapped_fragment(t,group) or any(bool(b.a in ids)!=bool(b.b in ids) for b in s.bonds)
    if not changed or bond_orders(s)==bond_orders(t):raise FormalError('Empty/repeated-resonance loop has no net product formation')
    from .smiles import to_rdkit
    from rdkit import Chem
    def unmapped_inventory(v):
        mol=to_rdkit(v)
        for a in mol.GetAtoms():a.SetAtomMapNum(0)
        return Chem.MolToSmiles(mol,canonical=True,isomericSmiles=True,kekuleSmiles=True)
    # This is an additional no-net-change guard. All atom-mapped identity,
    # carrier/electron/stereo checks above remain mandatory.
    if unmapped_inventory(s)==unmapped_inventory(t):raise FormalError('Role exchange alone is not nonzero net chemical turnover')
    receipt=hashlib.sha256(json.dumps({'record':document.record_id,'revision':document.revision,'events':cycle.event_ids,'source':mapped_fragment(s,[a.atom_id for a in s.atoms]),'target':mapped_fragment(t,[a.atom_id for a in t.atoms]),'products':sorted(seen)},sort_keys=True).encode()).hexdigest()
    return {'status':'closed_one_turnover','compressible':True,'receipt':receipt,'product_atom_ids':sorted(seen),'reward':1,
            'chemical_evidence':'declared_not_independently_approved','chemical_reward':None,'reward_scope':'formal_once_only_receipt_not_chemical_correctness','further_turnovers_executed':0}
def reward_once(result,receipts):
    if result.get('status')!='closed_one_turnover' or result['receipt'] in receipts:return 0
    receipts.add(result['receipt']);return 1
