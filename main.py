from proteinModel import proteinModel
from comparator import Comparator

WT_models = proteinModel.quick_load(
    {
     '/home/sulcjo/Desktop/myomedin/mmpbsa/4G6F/WT/***/xmgrace/***' : ( (2,6,8), ('rmsd_backbone-vs-start.xvg', 'rmsf_all_atom.xvg', 'gyrate.xvg', 'minimal-periodic-distance.xvg'), ('rmsd', 'rmsf', 'rg', 'mpd') ),
     '/home/sulcjo/Desktop/myomedin/mmpbsa_gmxmmpbsa/WT/***/results/FINAL_RESULTS_MMPBSA.dat' : ( (2,6,8), ('gmxmmpbsa') ),
     '/home/sulcjo/Desktop/myomedin/sirah_umbrella/sirah_umbrella_wt/***/wham_results/***' : ( (2,6,8), ('profile_errors.xvg', 'histo.xvg', 'contacts_pulling.xvg'), ('umbrella_profile', 'umbrella_histogram', 'contacts')),
     '/home/sulcjo/Desktop/myomedin/sirah_umbrella/sirah_umbrella_wt/***/model_***/***/***/***.pdb{2,3': ( (2,6,8), ('a','b'), ('sequence') )
    }, handles = {'4G6F-MyoWT mod. ***' : ((2,6,8)) }
)

# '/path/***/***/***/***{2,2': ( (2,6,8), (1,2,3), ('dataset name') )




WT_compare = Comparator(proteinModels=WT_models)
WT_compare.setLineColor = 'blue'
WT_compare.setFontSizeLarge = 22
WT_compare.setFontSizeMedium = 18
WT_compare.setFigSize = (20, 15)
#WT_compare.plot_simple_dataset('rmsd')
#WT_compare.plot_simple_dataset('rg')
#WT_compare.plot_simple_dataset('rmsf')
#WT_compare.plot_mmpbsa('gmxmmpbsa_pb', title='PB')
#WT_compare.plot_mmpbsa('gmxmmpbsa_gb', title='GB')
#WT_compare.plot_umbrella(fit=True, stderror=True)
#WT_compare.plot_simple_dataset('contacts')
#WT_compare.show()