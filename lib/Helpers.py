from Morph import Morph

def make_morph():
    args = dict()
    args['genome'] = '3'
    args['src_model'] = '5'
    args['probanno'] = '15'
    args['protcomp'] = '6'
    args['genomews'] = '9145'
    args['src_modelws'] = '9145'
    args['probannows'] = '9145'
    args['protcompws'] = '9145'
    return Morph(args)

