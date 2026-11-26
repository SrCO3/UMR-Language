"""Donor reservations and finite ionic L/X/X2 accounting, separate from N."""
from collections import Counter
from ..chemistry import FormalError, Unsupported, atoms_by_id, bond_orders, bond_valences
from .elements import METAL_GROUP, METAL_OXIDATION


def checks(view,*,_metal_group=METAL_GROUP,_metal_oxidation=METAL_OXIDATION,
           _transition_metals=frozenset({'Cr','Ru','Pd','Zn','Hg'})):
    METAL_GROUP, METAL_OXIDATION = _metal_group, _metal_oxidation
    atoms=atoms_by_id(view);orders=bond_orders(view);centers={m.atom_id:m for m in view.metal_centers}
    mids={i for i,a in atoms.items() if a.element in METAL_GROUP}
    if set(centers)!=mids:raise FormalError('Every supported metal atom needs exactly one MetalCenter')
    ligands={l.ligand_id:l for l in view.ligands};membership=set();covered=set();reserved=Counter();used=set();sites=set()
    for l in view.ligands:
        members=set(l.atom_ids)
        if not members<=atoms.keys() or members&mids or not set(l.donor_atom_ids)<=members:
            raise FormalError('Invalid ligand atom/donor references')
        if len(members)!=len(l.atom_ids) or membership&members:raise Unsupported('Overlapping/bridging ligand definitions unsupported')
        membership |= members
        todo=[min(members)];seen=set()
        while todo:
            i=todo.pop()
            if i in seen:continue
            seen.add(i);todo += [j for edge in orders if i in edge for j in edge if j in members and j not in seen]
        if seen!=members:raise FormalError('Ligand must be a connected fragment')
        if l.ionic_class!='unknown' and l.ionic_charge!={'L':0,'X':-1,'X2':-2}[l.ionic_class]:
            raise FormalError('Ionic ligand class/charge mismatch')
    for binding in view.metal_bindings:
        m=binding.metal_atom_id
        if m not in centers or binding.ligand_id not in ligands:raise FormalError('Missing metal/ligand binding reference')
        l=ligands[binding.ligand_id];members=set(l.atom_ids);donors=binding.donor_atom_ids
        if binding.ligand_id in used:raise Unsupported('Chelated/multimetal/bridging ligand scope not registered')
        used.add(binding.ligand_id)
        if len(set(donors))!=len(donors) or set(donors)!=set(l.donor_atom_ids):raise FormalError('Binding/ligand donor mismatch')
        if binding.site_id:
            if binding.site_id not in centers[m].site_ids or (m,binding.site_id) in sites:raise FormalError('Unknown or multiply occupied metal site')
            sites.add((m,binding.site_id))
        valence_to_metal=sum(orders.get(tuple(sorted((m,i))),0) for i in members)
        ionic_q=sum(atoms[i].formal_charge for i in members)-valence_to_metal
        if l.ionic_charge is not None and ionic_q!=l.ionic_charge:raise FormalError('Ligand charge differs from declared ionic partition')
        if binding.mode=='sigma_covalent':
            if len(donors)!=1 or binding.hapticity!=1:raise Unsupported('Only monodentate sigma X/X2 ligands implemented')
            pair=tuple(sorted((m,donors[0])));o=orders.get(pair,0)
            if o not in {1,2}:raise Unsupported('Only sigma single or formal terminal oxo double bonds supported')
            if l.ionic_class=='X2' and (atoms[donors[0]].element!='O' or len(members)!=1):raise Unsupported('X2 limited to terminal oxo')
            if l.ionic_class not in {'X','X2','unknown'}:raise Unsupported('Covalent L representation requires explicit conversion')
            covered.add(pair)
        elif binding.mode=='dative':
            if len(donors)!=1 or binding.hapticity!=1 or l.ionic_class not in {'L','unknown'}:raise Unsupported('Only monodentate two-electron L dative ligands')
            if valence_to_metal:raise FormalError('Dative pair counted also as covalent bond')
            reserved[donors[0]]+=2
        else:
            if binding.hapticity!=2 or len(donors)!=2 or l.ionic_class!='L' or any(atoms[i].element!='C' for i in donors) or orders.get(tuple(sorted(donors)))!=2:
                raise Unsupported('Eta scope is one neutral alkene pair')
            if valence_to_metal:raise FormalError('Eta unit duplicated as ordinary bonds')
        for i in members:
            if any(j not in members and j!=m for edge in orders if i in edge for j in edge if j!=i):
                raise FormalError('Bound ligand omits a covalently attached atom')
    metal_edges={p for p in orders if set(p)&mids}
    if any(set(p)<=mids for p in metal_edges):raise Unsupported('Metal-metal bonds require a dedicated profile')
    if covered!=metal_edges:raise FormalError('Every metal covalent bond requires a ligand binding')
    for adduct in view.lewis_adducts:
        d,a=adduct.donor_atom_id,adduct.acceptor_atom_id
        if atoms[a].element!='B' or d in mids or atoms[d].element not in {'N','O','P','S'}:
            raise Unsupported('Lewis-adduct scope: one N/O/P/S donor to trivalent B')
        if (bond_valences(view)[a],atoms[a].nonbonding_electrons,atoms[a].formal_charge)!=(3,0,0):
            raise FormalError('Lewis adduct requires neutral trivalent B with a vacant acceptor orbital')
        if orders.get(tuple(sorted((d,a)))):raise FormalError('Lewis adduct duplicated as covalent bond')
        if sum(x.acceptor_atom_id==a for x in view.lewis_adducts)>1:raise FormalError('B acceptor capacity exceeded')
        reserved[d]+=2
    for i,n in reserved.items():
        atom=atoms[i]
        if atom.radical_electrons is None or atom.nonbonding_electrons is None:raise Unsupported('Unknown donor pair occupation')
        if n>atom.nonbonding_electrons-atom.radical_electrons:raise FormalError('Repeated consumption of donor-owned coordination pair')
    counts=[]
    for i,m in centers.items():
        bs=[b for b in view.metal_bindings if b.metal_atom_id==i]
        members={i}|{j for b in bs for j in ligands[b.ligand_id].atom_ids}
        q=sum(atoms[j].formal_charge for j in members)
        known=all(ligands[b.ligand_id].ionic_charge is not None for b in bs)
        os=q-sum(ligands[b.ligand_id].ionic_charge for b in bs) if known else None
        el=atoms[i].element;g=METAL_GROUP[el]
        if os is not None and os not in METAL_OXIDATION[el]:raise Unsupported(f'{el} oxidation state {os} outside registered scope')
        d=(min(10,g-os) if el in _transition_metals else None) if os is not None else None
        ve=g-os+sum(4 if ligands[b.ligand_id].ionic_class=='X2' else 2 for b in bs) if known else None
        for field,value in [('oxidation_state',os),('d_electrons',d),('valence_electrons',ve)]:
            stored=getattr(m,field)
            if stored is not None and stored!=value:raise FormalError(f'Metal {i} {field} declaration differs from ionic analysis')
        if m.spin_state=='closed_shell' and atoms[i].radical_electrons not in {0,None}:raise FormalError('Metal closed-shell declaration contradicts r')
        if m.spin_state=='open_shell' and atoms[i].radical_electrons==0:raise FormalError('Metal open-shell declaration contradicts r')
        if len(set(m.site_ids))!=len(m.site_ids):raise FormalError('Duplicate coordination sites')
        if m.vacant_site_ids is not None and (not set(m.vacant_site_ids)<=set(m.site_ids) or any((i,s) in sites for s in m.vacant_site_ids)):
            raise FormalError('Vacant site is occupied or undeclared')
        own={b.binding_id for b in bs}
        if any(r.first_binding_id not in own or r.second_binding_id not in own or r.first_binding_id==r.second_binding_id for r in m.relative_positions):raise FormalError('Invalid relative-position reference')
        counts.append({'metal_atom_id':i,'element':el,'complex_atom_ids':sorted(members),'complex_total_charge':q,
                       'atomic_formal_charge':atoms[i].formal_charge,'oxidation_state':os,'d_electrons':d,
                       'valence_electrons':ve,'coordination_units':len(bs),'coordination_number':sum(len(b.donor_atom_ids) for b in bs),
                       'counting_method':'ionic','global_ledger_electrons_added':0,'eighteen_electron_rule':'advisory_only',
                       'spin_status':'unknown' if atoms[i].radical_electrons is None else 'declared_local_occupation',
                       'X2_scope':'terminal oxo as four-electron ionic donor; LX2 alternative not inferred'})
    return reserved,counts
