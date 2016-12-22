from perlfunc.perlfunc import perlfunc, perl5lib, perlreq, perlargs


@perlfunc
@perlreq("myfbaclient.pm")
def call(function, parameters):
    pass


def gapfill_model(args):
    info = call("gapfill_metabolic_model", args)
    return dict(zip(info[0::2], info[1::2]))


def translate_model(args):
    info = call("propagate_model_to_new_genome", args)
    return dict(zip(info[0::2], info[1::2]))


def add_reactions(args):
    info = call("edit_metabolic_model", args)
    return dict(zip(info[0::2], info[1::2]))


def edit_reactions(args):
    info = call("edit_metabolic_model", args)
    return dict(zip(info[0::2], info[1::2]))


def remove_reactions(args):
    call("edit_metabolic_model", args)


def compare_models(args):
    info = call("compare_models", args)
    return dict(zip(info[0::2], info[1::2]))


def build_metabolic_model(args):
    info = call("build_metabolic_model", args)
    return dict(zip(info[0::2], info[1::2]))


def runfba(args):
    info = call("run_flux_balance_analysis", args)
    return dict(zip(info[0::2], info[1::2]))

