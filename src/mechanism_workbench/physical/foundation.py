"""Structural preconditions used by the polar event checker."""
from ..umr.chemistry import FormalError, Unsupported, order, bond_valences


def check_input(template_id,view,roles):
    """Finite structural prerequisites; no claim about rates or pathway dominance."""
    atoms={a.atom_id:a for a in view.atoms};vals=bond_valences(view)
    def bonded(a,b,n=1):
        if order(view,roles[a],roles[b])!=n:raise FormalError(f'Template needs {a}-{b} bond order {n}')
    if template_id in {'SN1','SN2'}:
        c=roles['carbon'];a=atoms[c];bonded('carbon','leaving')
        if (vals[c],a.formal_charge,a.nonbonding_electrons)!=(4,0,0) or any(b.order!=1 for b in view.bonds if c in {b.a,b.b}):raise FormalError('Alkyl substitution requires saturated tetravalent carbon')
        if atoms[roles['leaving']].element not in {'F','Cl','Br','I','O','S'}:raise Unsupported('Unregistered alkyl leaving-group class')
        if template_id=='SN2' and sum(atoms[b.b if b.a==c else b.a].element=='C' for b in view.bonds if c in {b.a,b.b})>2:raise Unsupported('This SN2 template excludes tertiary carbon')
    if template_id in {'CARBONYL_ADDITION','ACYL_SUBSTITUTION','TAUTOMER_BASE','TAUTOMER_ACID'}:
        bonded('carbon','oxygen',2)
    if template_id=='ACYL_SUBSTITUTION':
        bonded('carbon','leaving')
        if atoms[roles['leaving']].element in {'C','H'}:raise Unsupported('Acyl template does not eject carbon or hydride')
    if template_id.startswith('TAUTOMER_'):
        bonded('alpha','carbon');bonded('alpha','hydrogen')
    if template_id in {'SNAR_STEPWISE','SNAR_CONCERTED','EAS'}:
        start=roles['ipso'];cycles=[]
        def visit(path):
            if len(path)==6:
                if order(view,path[-1],start)==1:cycles.append(path)
                return
            for b in view.bonds:
                if path[-1] not in {b.a,b.b}:continue
                j=b.b if b.a==path[-1] else b.a
                if j not in path and atoms[j].element=='C' and b.order==(2 if len(path)%2 else 1):visit(path+[j])
        visit([start])
        if not cycles or ('adjacent' in roles and not any(p[1]==roles['adjacent'] for p in cycles)):raise Unsupported('This aromatic template needs a six-carbon alternating Lewis ring and declared pi partner')
        if template_id.startswith('SNAR'):
            bonded('ipso','leaving')
            if atoms[roles['leaving']].element not in {'F','Cl','Br','I'}:raise Unsupported('Initial SNAr template restricted to aryl halides')
        else:bonded('ipso','hydrogen')
