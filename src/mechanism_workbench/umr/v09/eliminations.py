"""Distinct finite concerted topologies, not a universal multiple-bond permission."""
from ..ast import Shift
from ..chemistry import atoms_by_id,order,FormalError,Unsupported
SPECS={'XANTHATE_SYN_ELIMINATION':['sulfur','carbon','oxygen','alpha','beta','hydrogen'],
       'CARBON_ELECTROCYCLIC_6PI_OPEN':['a','b','c','d','e','f'],
       'OXAPHOSPHETANE_CYCLOREVERSION':['a','b','phosphorus','oxygen']}
def compile(v,e):
    r=e.general_roles;ids=[r[k] for k in SPECS[e.general_event]];ats=atoms_by_id(v)
    def bond(a,b,n):
        if order(v,a,b)!=n:raise FormalError(f'Expected {a}-{b} order {n}')
    if e.general_event=='CARBON_ELECTROCYCLIC_6PI_OPEN':
        a,b,c,d,f,g=ids
        if any(ats[i].element!='C' or ats[i].formal_charge or ats[i].radical_electrons for i in ids):
            raise Unsupported('Finite neutral closed-shell six-carbon electrocyclic opening')
        for x,y,n in [(a,b,1),(b,c,2),(c,d,1),(d,f,2),(f,g,1),(g,a,1)]:bond(x,y,n)
        if e.activation!='thermal':raise Unsupported('This opening requires an explicit thermal channel')
        if e.rotation_mode not in {'unknown','disrotatory'}:raise FormalError('Thermal six-electron opening requires disrotatory compatibility')
        return [Shift(a=g,b=a,c=b),Shift(a=b,b=c,c=d),Shift(a=d,b=f,c=g)],[]
    if e.general_event=='XANTHATE_SYN_ELIMINATION':
        s,c,o,a,b,h=ids
        if [ats[i].element for i in ids]!=['S','C','O','C','C','H']:raise Unsupported('Xanthate S=C-O-C-C-H topology required')
        for x,y,n in [(s,c,2),(c,o,1),(o,a,1),(a,b,1),(b,h,1),(h,s,0)]:bond(x,y,n)
        if e.activation not in {'thermal','unknown'}:raise Unsupported('Finite thermal syn Ei channel; photochemistry distinct')
        return [Shift(a=c,b=s,c=h),Shift(a=h,b=b,c=a),Shift(a=a,b=o,c=c)],[]
    a,b,p,o=ids
    if [ats[i].element for i in ids]!=['C','C','P','O']:raise Unsupported('Four-center C-C-P-O cycloreversion')
    for x,y in [(a,b),(a,p),(p,o),(o,b)]:bond(x,y,1)
    return [Shift(a=p,b=a,c=b),Shift(a=b,b=o,c=p)],[]
