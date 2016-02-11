
#import WCM_tools

def TRSL_module_dict():
    """
    Computational model of S. cerevisiae cell cycle.
    TRSL: nonspecific discrete protein translation model.
    Python module written by Martin Seeger.
    """

    info = "Computational model of S. cerevisiae cell cycle. TRSL: nonspecific discrete protein translation model. Python module written by Martin Seeger."
    
    module_dict = {}

    ### Model Name:
    module_dict['name'] = 'TRSL: discrete protein translation'

    ### Initial Values for Species:
    module_dict['initvars'] = {"protein": 0, "ribos._bound": 0, "ribos._free": 200000, "tRNA_bound": 0,
                               "tRNA_free": None, 'GTP': None, 'GDP': None, 'ATP': None, 'AMP': None}  #TODO: FIXME: values

    ### Initial Values for Parameters:
    module_dict['initpars'] = {'Kp_Sic1': 15.028, 'n1': 6.0}  #TODO: FIXME: parameters

    ### Model Parameters:
    module_dict['pars'] = ['Kp_Sic1', 'n1']  #TODO: FIXME: parameters

    ### Model Species:
    module_dict['vars'] = ["protein", "ribos._bound", "ribos._free", "tRNA_bound", "tRNA_free", "ATP", "AMP", "GTP", "GDP"]

    ### Species Units:
    module_dict['units'] = {'protein': 'dimensionless', 'ribos._bound': 'dimensionless', 'ribos._free': 'dimensionless',
                            'tRNA_bound':'dimensionless', 'GTP': 'dimensionless', 'GDP': 'dimensionless',
                            'ATP': 'dimensionless', 'AMP': 'dimensionless'}  #TODO: FIXME: units

    ### Species Compartment:
    module_dict['sp_compartment'] = {"protein": 'cytosol', "ribos._bound": 'cytosol', "ribos._free": 'cytosol',
                                     "tRNA_bound": 'cytosol', "tRNA_free": 'cytosol', 'GTP': 'cytosol',
                                     'GDP': 'cytosol', 'ATP': 'cytosol', 'AMP': 'cytosol'}

    ### Compartment Annotations:
    module_dict['com_annotations'] = {'cytosol': 'GO:0005829',
                                      #'nucleus': 'GO:0005634'
                                      }

    ### Species Annotations:
    module_dict['sp_annotations'] = {"protein": "CHEBI:36080",  # generic protein
                                     # "amino_acid": "CHEBI:15705", # generic amino acid
                                     "ribos._bound": "GO:0042788",
                                     "ribos._free": "GO:0005840",
                                     "tRNA_bound": "CHEBI:17843_b",
                                     "tRNA_free": "GO:0005564",  # generic tRNA; http://www.ebi.ac.uk/chebi/searchId.do;010D9AC7FDDC72158F86B943C40AD04A?chebiId=CHEBI:2651 lists some others
                                     'GTP': 'CHEBI:15996', #TODO: FIXME: are GTP and ATP different in the models?
                                     'GDP': 'CHEBI:17552',
                                     'ATP': 'CHEBI:15422',
                                     'AMP': 'CHEBI:16027'}

    ### States:
    module_dict['sp_states'] = {var: '0' for var in module_dict['vars']}

    module_dict['units_pars'] = {'K_MBF': 'dimensionless', 'v0_Mcm1': 's**-1'}

    ### Reactions:
    reactions = {}
    reactions['v_MBF_akt'] = {'rate': {}, 'products': {}, 'substrates': {}, 'modifiers': {}}
    reactions['v_Cln3_p'] = {'rate': {}, 'products': {}, 'substrates': {}, 'modifiers': {}}

    ### Rates:
    # v_MBF_akt
    reactions['v_MBF_akt']['rate'] = '(     ( ( kp_MBF  *  ( Cln2_cyt  **  n1 ) )  /  ( ( K_MBF  **  n1 )  + ( Cln2_cyt  **  n1 ) ) )   )'
    # v_Cln3_p
    reactions['v_Cln3_p']['rate'] = '(    kp_Cln3  )'


    ### Substrates:
    # v_MBF_akt
    reactions['v_MBF_akt']['substrates'] = {}
    # v_Cln3_p
    reactions['v_Cln3_p']['substrates'] = {}


    ### Products:
    # v_MBF_akt
    reactions['v_MBF_akt']['products'] = {'MBF_nuc': 1.0}
    # v_Cln3_p
    reactions['v_Cln3_p']['products'] = {'Cln3_cyt': 1.0}


    ### Modifiers:
    # v_MBF_akt
    reactions['v_MBF_akt']['modifiers'] = {'Cln2_cyt': 1.0}
    # v_Cln3_p
    reactions['v_Cln3_p']['modifiers'] = {}

    module_dict['functions'] = {}
    module_dict['alg_eqs'] = {}
    module_dict['info'] = info
    module_dict['reactions'] = reactions
    #module_dict['odes'] = WCM_tools.generateODEsFromReactions(reactions, module_dict['sp_compartment'], module_dict['sp_annotations'], module_dict['units'])

    return module_dict
