from perlfunc.perlfunc import perlfunc, perl5lib, perlreq, perlargs


@perlfunc
@perlreq("myfbaclient.pm")
def call(function, parameters):
    pass


def gapfill_model(args):
    call("gapfill_metabolic_model", args)


def translate_model(args):
    call("propagate_model_to_new_genome", args)


def add_reactions(args):
    call("edit_metabolic_model", args)


def edit_reactions(args):
    call("edit_metabolic_model", args)


def remove_reactions(args):
    call("edit_metabolic_model", args)


def compare_models(args):
    call("compare_models", args)


def build_metabolic_model(args):
    call("build_metabolic_model", args)


def runfba(args):
    call("run_flux_balance_analysis", args)

