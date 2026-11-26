from typing import Literal
from pydantic import Field, model_validator
from ..ast import AtomID, Evidence, Identifier, Node
from ..v03.ast import Ligand as Ligand03, MetalCenter as MetalCenter03, MetalBinding
from ..v04.ast import (Atom, ChemicalState, ElectronicState, ReactionEvent as Event04,
                       LewisView as View04, MechanismDocument as Document04,
                       AnalysisRecord as Analysis04, ResonanceTransform)


class ElementDeclaration(Node):
    atom_id: AtomID
    valence_category: str = 'unknown'  # descriptive category, never a charge override
    oxidation_state: int | None = None
    coordination_number: int | None = Field(default=None, ge=0)
    evidence: list[Evidence] = Field(default_factory=list)


class Ligand(Ligand03):
    ionic_class: Literal['L', 'X', 'X2', 'unknown'] = 'unknown'
    assignment_basis: str = 'unknown'


class MetalCenter(MetalCenter03):
    spin_state: Literal['closed_shell', 'open_shell', 'unknown'] = 'unknown'


class LewisAdduct(Node):
    binding_id: Identifier
    donor_atom_id: AtomID
    acceptor_atom_id: AtomID
    electrons: Literal[2] = 2
    convention: Literal['donor_owned_dative/1'] = 'donor_owned_dative/1'
    evidence: list[Evidence] = Field(default_factory=list)


class LewisView(View04):
    element_rule_version: Literal['umr-elements/0.5'] = 'umr-elements/0.5'
    element_declarations: list[ElementDeclaration] = Field(default_factory=list)
    metal_centers: list[MetalCenter] = Field(default_factory=list)
    ligands: list[Ligand] = Field(default_factory=list)
    lewis_adducts: list[LewisAdduct] = Field(default_factory=list)

    @model_validator(mode='after')
    def element_references(self):
        ids={a.atom_id for a in self.atoms}
        declarations=[d.atom_id for d in self.element_declarations]
        if len(set(declarations))!=len(declarations) or not set(declarations)<=ids:
            raise ValueError('Element declarations require unique existing atoms')
        keys=[]
        for b in self.lewis_adducts:
            if b.donor_atom_id==b.acceptor_atom_id or not {b.donor_atom_id,b.acceptor_atom_id}<=ids:
                raise ValueError('Adduct endpoints must be distinct existing atoms')
            keys.append((b.donor_atom_id,b.acceptor_atom_id))
        if len(set(keys))!=len(keys) or len({b.binding_id for b in self.lewis_adducts})!=len(keys):
            raise ValueError('Duplicate adduct/electron-pair reservation')
        return self


class ReactionEvent(Event04):
    element_template: str | None = None


class AnalysisRecord(Analysis04):
    operation: Literal['RESONANCE_SCAN','SITE_ASSESSMENT','PKA_QUERY','REPRESENTATION_SENSITIVITY',
                       'METAL_ELECTRON_COUNT','METAL_REPRESENTATION_CHECK','ELEMENT_RULE_CHECK',
                       'ELEMENT_REPRESENTATION_CHECK','PHYSICAL_PARAMETER_QUERY']
    rule_version: str = 'umr/0.5'


class MechanismDocument(Document04):
    schema_version: Literal['umr/0.5'] = 'umr/0.5'
    views: list[LewisView] = Field(min_length=1)
    events: list[ReactionEvent] = Field(default_factory=list)
    analyses: list[AnalysisRecord] = Field(default_factory=list)
