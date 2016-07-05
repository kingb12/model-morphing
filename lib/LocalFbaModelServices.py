from perlfunc.perlfunc import perlfunc, perl5lib, perlreq, perlargs


@perlfunc
@perlreq("myfbaclient.pm")
def _call(function, parameters):
    pass


def gapfill_model(args):
    _call("gapfill_metabolic_model", args)


