from Chemical import *

class ChemicalReactionError(Exception):
    pass

class ChemicalReaction:

    def __init__(self, name:str, reactants:list[Chemical], products:list[Chemical], stoich_coeff:list[float], Q_reaction:float):
        '''
        Defines and holds information about chemical reactions. Performs reactions based on stoicheometry.
        '''
        if ((len(reactants) + len(products)) != len(stoich_coeff)): raise ChemicalReactionError("Number of chemical reaction coefficient must match number of reactants and products!")
        self.name = name
        self.reactants = reactants # array-like of valid coolprop chemical names
        self.products = products # array-like of valid coolprop chemical names
        self.reactant_stoich_coeff = stoich_coeff[:len(reactants)]
        self.product_stoich_coeff = stoich_coeff[len(reactants):]
        self.Q_reaction= Q_reaction

    def complete_reaction(self, reactant_flowrates:list[float]) -> tuple[list[float], list[float], float]:
        '''
        Ideal chemical reaction class, based on stoicheometry.

        Parameters:
            reactant_flowrates: List of reactant flow rates (kg/s)

        Returns:
            reactant_flowrates: List of the remaining reactants flow rates (kg/s) after the reaction
            product_flowrates: List of the products flow rates (kg/s) after the reaction
            P_total: total power of the chemical reaction (W)
        '''
        if (len(self.reactants) != len(reactant_flowrates)): raise ChemicalReactionError("Must specify a reactant flow rate for each reactant!")
        
        reactant_molflowrates = [0] * len(self.reactants)
        stoichmetric_molflowrates = [0] * len(self.reactants)
        for k in range(len(self.reactants)):
            reactant_molflowrates[k] = reactant_flowrates[k] / (self.reactants[k].get_mol_weight()/1000)
            stoichmetric_molflowrates[k] = reactant_molflowrates[k] / self.reactant_stoich_coeff[k]

        molflowrate = np.min(stoichmetric_molflowrates) 

        for k in range(len(self.reactants)):
            reactant_molflowrates[k] = reactant_molflowrates[k] - molflowrate * self.reactant_stoich_coeff[k]
            reactant_flowrates[k] = (reactant_molflowrates[k] * self.reactants[k].get_mol_weight()/1000)

        product_molflowrates = [0] * len(self.products)
        product_flowrates = [0] * len(self.products)
        for k in range(len(self.products)):
            product_molflowrates[k] = molflowrate * self.product_stoich_coeff[k]
            product_flowrates[k] = (product_molflowrates[k] * self.products[k].get_mol_weight()/1000)

        P_total = (self.Q_reaction * molflowrate)
        return reactant_flowrates, product_flowrates, P_total
    
class ChemicalReactionSet:

    def __init__(self):
        self.reactions: dict[str, ChemicalReaction] = {}

    def set_reaction(self, aChemicalReaction:ChemicalReaction):
        self.reactions[aChemicalReaction.name] = aChemicalReaction

    def append_chemical_reaction_list(self, aChemicalReactionList:list[ChemicalReaction]): 
        # append multiple new chemical reactions to the set
        for k in range(len(aChemicalReactionList)):
            self.reactions[aChemicalReactionList[k].name] = aChemicalReactionList[k] 

    def get_reaction(self, reaction_name:str) -> ChemicalReaction:
        return self.reactions[reaction_name]