from perlfunc.perlfunc import perlfunc, perl5lib, perlreq, perlargs


@perlfunc
@perlreq("myfbaclient.pm")
def call(function, parameters):
    pass


def gapfill_model(args):
    return zip(call("gapfill_metabolic_model", args))


def translate_model(args):
    return zip(call("propagate_model_to_new_genome", args))


def add_reactions(args):
    return zip(call("edit_metabolic_model", args))


def edit_reactions(args):
    return zip(call("edit_metabolic_model", args))


def remove_reactions(args):
    return zip(call("edit_metabolic_model", args))


def compare_models(args):
    return zip(call("compare_models", args))


def build_metabolic_model(args):
    return zip(call("build_metabolic_model", args))


def runfba(args):
    return zip(call("run_flux_balance_analysis", args))

