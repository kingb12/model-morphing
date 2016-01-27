from Helpers import *
m = make_morph()
m = Client.prepare_supermodel(m, fill_src=True)
print runfba(m)
rxns = Client.removal_list(m.rxn_labels['gene-match'])
m = Client.process_reactions(m, rxns)
dump(m, '../data/morph-gene-match-removal.pkl')
