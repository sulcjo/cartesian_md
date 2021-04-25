from proteinModel import proteinModel
from comparator import Comparator

base_path = '/home/sulcjo/Desktop/myomedin'
WT_models = proteinModel.quick_load(
    {
     f'{base_path}/mmpbsa/4G6F/WT/***/xmgrace/***' : ( (2,6,8), ('rmsd_backbone-vs-start.xvg', 'rmsf_all_atom.xvg', 'gyrate.xvg', 'minimal-periodic-distance.xvg'), ('rmsd', 'rmsf', 'rg', 'mpd') ),
     f'{base_path}/mmpbsa_gmxmmpbsa/WT/***/results/FINAL_RESULTS_MMPBSA.dat' : ( (2,6,8), ('gmxmmpbsa') ),
     f'{base_path}/sirah_umbrella/sirah_umbrella_wt/***/wham_results/***' : ( (2,6,8), ('profile_errors.xvg', 'histo.xvg', 'contacts_pulling.xvg'), ('umbrella_profile', 'umbrella_histogram', 'contacts')),
     f'{base_path}/sirah_umbrella/sirah_umbrella_wt/***/model_***.pdb^2': ( (2,6,8), ('sequence') )
    }, handles = {'4G6F-MyoWT mod. ***' : ((2,6,8)) }
)

# '/path/***/***/***/***{2,2': ( (2,6,8), (1,2,3), ('dataset name') )

WT_models[1].split = [(0,239),(240,350)]

WT_models[1].get_GRINN_datasets('/home/sulcjo/Desktop/myomedin/grinn_linux_v110_hf1/WT_10E8_model_6/')

#print(WT_models[1].pdb_file)


#WT_models[1].DSSP_assign()
#WT_models[1].call_pymol()


WT_compare = Comparator(proteinModels=WT_models)
WT_compare.setLineColor = 'blue'
WT_compare.setFontSizeLarge = 22
WT_compare.setFontSizeMedium = 18
WT_compare.setFigSize = (20, 15)
WT_compare.plot_simple_dataset('rmsd')
WT_compare.plot_simple_dataset('rg')
WT_compare.plot_simple_dataset('rmsf')
WT_compare.plot_mmpbsa('gmxmmpbsa_pb', title='PB')
WT_compare.plot_mmpbsa('gmxmmpbsa_gb', title='GB')
WT_compare.plot_umbrella(fit=True, stderror=True)
WT_compare.plot_simple_dataset('contacts')
#WT_compare.show()