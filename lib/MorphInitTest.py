from Morph import Morph
a = Morph()
for i in dir(a):
    print i
    print getattr(a, i)
