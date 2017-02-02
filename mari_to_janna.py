import sys
import model_morphing.Helpers as Helpers
import model_morphing.objects as obj
morph = Helpers.mari_to_janna(int(sys.argv[1]))
morph.prepare_supermodel(fill_src=False)
new_media = obj.Media(25, 9145)
morph.model = obj.FBAModel(6, 14620)
morph.translate_media(new_media)
morph.process_reactions()
Helpers.dump(morph, '../data/mari_to_janna_morph_revised.pkl')
