import sys
import model_morphing.Helpers as Helpers
from model_morphing.objects import Media, FBAModel
import model_morphing.GrowthConditions as GrowthConditions
#runs.mari_to_bark_morph(13227, filename='../data/mari_to_bark_h2co2_morph.pkl')
#runs.mari_to_bark_morph(13228, newmedia=7, newmediaws=12091,filename='../data/mari_to_bark_meoh_morph.pkl')
#runs.mari_to_bark_morph(13229, newmedia=8, newmediaws=12091,filename='../data/mari_to_bark_ac_morph.pkl')
morph = Helpers.mari_to_bark(int(sys.argv[1]))
morph.prepare_supermodel(fill_src=False)
media = [Media(24, 9145), Media(7, 12091), Media(8, 12091)]
for m in media:
    morph.translate_media(m)
morph.process_reactions(growth_condition=GrowthConditions.AllMedia(media))
Helpers.dump(morph, '../data/mari_to_bark_3media_morph.pkl')

