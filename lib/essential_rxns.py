from biokbase.fbaModelServices.Client import fbaModelServices
reactions = ['rxn05380', 'rxn05364', 'rxn05397', 'rxn05401', 'rxn05393', 'rxn05385', 'rxn05368', 'rxn05360', 'rxn05372', 'rxn05376', 'rxn05405', 'rxn05389']
with open (".kbase_fbaModelServicesURL", "r") as myfile:
    url = myfile.read().replace('\n','')
fbac = fbaModelServices(url)
rxns = fbac.get_reactions({'reactions' : reactions})
print rxns[0].keys()
for r in rxns:
    print r['definition']
