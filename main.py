from proteinModel import proteinModel
from comparator import Comparator

base_path = '/home/sulcjo/Desktop/myomedin'
"""
all = proteinModel.quick_load(
    {
     f'{base_path}/mmpbsa/4G6F/***/***/xmgrace/***' : ( (24,25,92,158,'WT','MPER'),(range(0,11)), ('rmsd_backbone-vs-start.xvg', 'rmsf_all_atom.xvg', 'gyrate.xvg', 'minimal-periodic-distance.xvg'), ('rmsd', 'rmsf', 'rg', 'mpd') ),
     f'{base_path}/mmpbsa_gmxmmpbsa/***/***/results/FINAL_RESULTS_MMPBSA.dat' : ( (24,25,92,158,'WT','MPER'),(range(0,11)), ('gmxmmpbsa') ),
    }, handles = {'4G6F-Myo*** mod. ***' : ((24,25,92,158,'WT','MPER'),(range(0,11)))}
)
all_compare = Comparator(proteinModels=all)
all_compare.setLineColor = 'blue'
all_compare.setFontSizeLarge = 22
all_compare.setFontSizeMedium = 18
all_compare.setFigSize = (20, 15)
#all_compare.plot_simple_dataset('rmsd')
#all_compare.plot_simple_dataset('rg')
all_compare.plot_mmpbsa('gmxmmpbsa_pb',title='PB')
all_compare.show()
#WT_compare.plot_simple_dataset('rg')
#WT_compare.plot_simple_dataset('rmsf')
#WT_compare.plot_mmpbsa('gmxmmpbsa_pb', title='PB')

"""
"""
myo24 = proteinModel.quick_load(
    {
     f'{base_path}/mmpbsa/4G6F/24/***/xmgrace/***' : ( (0,4), ('rmsd_backbone-vs-start.xvg', 'rmsf_all_atom.xvg', 'gyrate.xvg', 'minimal-periodic-distance.xvg'), ('rmsd', 'rmsf', 'rg', 'mpd') ),
     f'{base_path}/mmpbsa_gmxmmpbsa/WT/***/results/FINAL_RESULTS_MMPBSA.dat' : ( (0,4), ('gmxmmpbsa') ),
     f'{base_path}/sirah_umbrella/sirah_umbrella_24/***/wham_results/***' : ( (0,4), ('profile_errors.xvg', 'histo.xvg', 'contacts_pulling.xvg'), ('umbrella_profile', 'umbrella_histogram', 'contacts')),
     f'{base_path}/mmpbsa_gmxmmpbsa/24/***/grinn_output/system_dry.pdb': ( (0,4), ('sequence') )
    }, handles = {'4G6F-Myo24 mod. ***' : ((0,4)) }
)

# '/path/***/***/***/***{2,2': ( (2,6,8), (1,2,3), ('dataset name') )

myo24[0].split = [(0,239),(240,350)]
myo24[1].split = [(0,239),(240,350)]
myo24[0].DSSP_assign()
myo24[1].DSSP_assign()
myo24[1].get_GRINN_datasets('/home/sulcjo/Desktop/myomedin/mmpbsa_gmxmmpbsa/24/4/grinn_output/')
myo24[0].get_GRINN_datasets('/home/sulcjo/Desktop/myomedin/mmpbsa_gmxmmpbsa/24/0/grinn_output/')
#WT_models[1].call_pymol()
#print(WT_models[1].pdb_file)


#WT_models[1].DSSP_assign()
#WT_models[1].call_pymol()


myo24_compare = Comparator(proteinModels=myo24)
myo24_compare.setLineColor = 'blue'
myo24_compare.setFontSizeLarge = 22
myo24_compare.setFontSizeMedium = 18
myo24_compare.setFigSize = (20, 15)
#WT_compare.plot_simple_dataset('rmsd')
#WT_compare.plot_simple_dataset('rg')
#WT_compare.plot_simple_dataset('rmsf')
#WT_compare.plot_mmpbsa('gmxmmpbsa_pb', title='PB')
#WT_compare.plot_mmpbsa('gmxmmpbsa_gb', title='GB')
#WT_compare.plot_umbrella(fit=True, stderror=True, check_sampling=True)
#WT_compare.plot_simple_dataset('contacts')
mark_resis=[pos+(240-14)-1 for pos in (34,35,36,37,63,64,65,89,90,91,92,93)]
#WT_compare.plot_IEM(modelIndexes=[1], dataset='vdw_IEM', mark_resis=mark_resis)
#WT_compare.plot_IEM(modelIndexes=[0], dataset='total_IEM', mark_resis=mark_resis, )
model01_delta_dataframe = myo24_compare.compare_IEM(modelIndexes=[0,1])
print(model01_delta_dataframe.max())

myo24_compare.show()

#mod24_2 = proteinModel(annotation='24_2')
#mod24_2.get_dataset(path='/home/sulcjo/Desktop/myomedin/sirah_umbrella/sirah_umbrella_24/2/wham_results/profile_errors.xvg', dataset_name='umbrella_profile')
#mod24_2.get_dataset(path='/home/sulcjo/Desktop/myomedin/sirah_umbrella/sirah_umbrella_24/2/wham_results/histo.xvg',dataset_name='umbrella_histogram')
#mod24_2_comp = Comparator(proteinModels=[mod24_2])
#mod24_2_comp.plot_umbrella(stderror=False, check_sampling=True, check_sampling_limit=0)
#mod24_2_comp.show()


sirah_umbrella = proteinModel.quick_load(
    {
        f'{base_path}/sirah_umbrella/sirah_umbrella_mper/***_unrestrained/wham_results/***' : ((2,6,'original'), ('profile_errors.xvg','histo.xvg','contacts_pulling.xvg'), ('umbrella_profile','umbrella_histogram','contacts')),

    }, handles= {'10E8-MPER mod. ***' : ((2,6,'original'))}

)

compare_sirah = Comparator(proteinModels=sirah_umbrella)
compare_sirah.plot_umbrella(stderror=True,fit=True,check_sampling=False)
compare_sirah.plot_simple_dataset(dataset='contacts')
compare_sirah.show()

atomistic_umbrella = proteinModel.quick_load(
    {
        f'{base_path}/umbrella_sampling/round_2_out/158_***/wham_results/***' : ((2,3,8), ('profile_errors.xvg','histo.xvg','contacts_pulling.xvg'), ('umbrella_profile','umbrella_histogram','contacts')),

    }, handles= {'10E8-158 mod. ***' : ((2,3,8))}

)

compare_atomistic = Comparator(proteinModels=atomistic_umbrella)
compare_atomistic.plot_umbrella(stderror=True,fit=False,max_traj=5,check_sampling=False)
compare_atomistic.plot_simple_dataset(dataset='contacts')
compare_atomistic.show()


sirah_umbrella = proteinModel.quick_load(
    {
        f'{base_path}/sirah_umbrella/sirah_umbrella_***/***_unrestrained/wham_results/***' : ((24,25,92,158,'wt','mper'),(0,1,2,3,4,5,6,7,8,9,'original'), ('profile_errors.xvg','histo.xvg','contacts_pulling.xvg'), ('umbrella_profile','umbrella_histogram','contacts')),

    }, handles= {'10E8-*** mod. ***' : ((24,25,92,158,'wt','mper'),(0,1,2,3,4,5,6,7,8,9,'original'))}

)

compare_sirah = Comparator(proteinModels=sirah_umbrella)
compare_sirah.plot_umbrella(stderror=True,fit=False,check_sampling=False, max_traj=5)
compare_sirah.plot_simple_dataset(dataset='contacts')
compare_sirah.show()

MPERs = proteinModel.quick_load(
    {
    f'{base_path}/mmpbsa_gmxmmpbsa/MPER/***/grinn_output/system_dry.pdb' : ((2,6,'original'), ('sequence')),
    }, handles={'10E8-MPER mod. ***' : (2,6,'original')}
)

#mark_resis=[pos+(240-14)-1 for pos in (34,35,36,37,63,64,65,89,90,91,92,93)]
for i,j in zip([0,1,2],[2,6,'original']):
    MPERs[i].split = [(0,238),(239,266)]
    #myo158[i].DSSP_assign()
    MPERs[i].get_GRINN_datasets(f'{base_path}/mmpbsa_gmxmmpbsa/MPER/{j}/grinn_output')

compareMPERs = Comparator(proteinModels=MPERs)
#model23_delta = compare158.compare_IEM(modelIndexes=[0,1])
compareMPERs.plot_IEM(modelIndexes=[0,1,2])
orig_6 = compareMPERs.compare_IEM(modelIndexes=[1,2])
compareMPERs.show()

MPER = proteinModel.quick_load(
    {
    f'{base_path}/mmpbsa_gmxmmpbsa/MPER/***/results/FINAL_RESULTS_MMPBSA.dat' : ((2,6,'original'), ('gmxmmpbsa')),
    f'{base_path}/mmpbsa_gmxmmpbsa/MPER/***/grinn_output/system_dry.pdb' : ((2,6,'original'), ('sequence')),
    }, handles={'10E8-MPER mod. ***' : (2,6,'original')}
)

#mark_resis=[pos+(240-14)-1 for pos in (34,35,36,37,63,64,65,89,90,91,92,93)]
MPER[0].split = [(0,130),(239,267)]
MPER[0].DSSP_assign()
MPER[0].get_GRINN_datasets(f'{base_path}/mmpbsa_gmxmmpbsa/MPER/2/grinn_output')
MPER[1].split = [(0,130),(239,267)]
MPER[1].DSSP_assign()
MPER[1].get_GRINN_datasets(f'{base_path}/mmpbsa_gmxmmpbsa/MPER/6/grinn_output')
MPER[2].get_dataset('/home/sulcjo/Desktop/myomedin/docking/cluspro_in/4G6F_original_mper.pdb','sequence')
MPER[2].split = [(0,130),(239,267)]
MPER[2].DSSP_assign()
MPER[2].get_GRINN_datasets(f'{base_path}/mmpbsa_gmxmmpbsa/MPER/original/grinn_output')



compareMPER = Comparator(proteinModels=MPER)
#model23_delta = compare158.compare_IEM(modelIndexes=[0,1])
compareMPER.plot_IEM(modelIndexes=[0,1,2], write_best_pairs=False, dataset='elec_IEM')
compareMPER.show()


"""





MLA024 = proteinModel.quick_load(
    {
    f'{base_path}/mmpbsa_gmxmmpbsa/24/***/results/FINAL_RESULTS_MMPBSA.dat' : ((2,10), ('gmxmmpbsa')),
    f'{base_path}/mmpbsa_gmxmmpbsa/24/***/grinn_output/system_dry.pdb' : ((2,10), ('sequence')),
    }, handles={'10E8-MLA024 mod. ***' : (2,10)}
)

mark_resis=[pos+(240-14)-1 for pos in (34,35,36,37,63,64,65,89,90,91,92,93)]
MLA024[0].split = [(0,130),(240,350)]
MLA024[0].DSSP_assign()
MLA024[0].get_GRINN_datasets(f'{base_path}/mmpbsa_gmxmmpbsa/24/2/grinn_output')


compareMLA024 = Comparator(proteinModels=MLA024)
#model23_delta = compare158.compare_IEM(modelIndexes=[0,1])
compareMLA024.plot_IEM(modelIndexes=[0], write_best_pairs=False, mark_resis=mark_resis, dataset='total_IEM')
compareMLA024.show()



MLA025 = proteinModel.quick_load(
    {
    f'{base_path}/mmpbsa_gmxmmpbsa/25/***/results/FINAL_RESULTS_MMPBSA.dat' : ((5,10), ('gmxmmpbsa')),
    f'{base_path}/mmpbsa_gmxmmpbsa/25/***/grinn_output/system_dry.pdb' : ((5,10), ('sequence')),
    }, handles={'10E8-MLA025 mod. ***' : (5,10)}
)

mark_resis=[pos+(240-14)-1 for pos in (34,35,36,37,63,64,65,89,90,91,92,93)]
MLA025[0].split = [(0,130),(240,350)]
MLA025[0].DSSP_assign()
MLA025[0].get_GRINN_datasets(f'{base_path}/mmpbsa_gmxmmpbsa/25/5/grinn_output')


compareMLA025 = Comparator(proteinModels=MLA025)
#model23_delta = compare158.compare_IEM(modelIndexes=[0,1])
compareMLA025.plot_IEM(modelIndexes=[0], write_best_pairs=False, mark_resis=mark_resis, dataset='total_IEM')
compareMLA025.show()


MLA092 = proteinModel.quick_load(
    {
    f'{base_path}/mmpbsa_gmxmmpbsa/92/***/results/FINAL_RESULTS_MMPBSA.dat' : ((0,10), ('gmxmmpbsa')),
    f'{base_path}/mmpbsa_gmxmmpbsa/92/***/grinn_output/system_dry.pdb' : ((0,10), ('sequence')),
    }, handles={'10E8-MLA092 mod. ***' : (0,10)}
)

mark_resis=[pos+(240-14)-1 for pos in (34,35,36,37,63,64,65,89,90,91,92,93)]
MLA092[0].split = [(0,130),(240,350)]
MLA092[0].DSSP_assign()
MLA092[0].get_GRINN_datasets(f'{base_path}/mmpbsa_gmxmmpbsa/92/0/grinn_output')

compareMLA092 = Comparator(proteinModels=MLA092)
#model23_delta = compare158.compare_IEM(modelIndexes=[0,1])
compareMLA092.plot_IEM(modelIndexes=[0], write_best_pairs=False, mark_resis=mark_resis, dataset='total_IEM')
compareMLA092.show()



MLA158 = proteinModel.quick_load(
    {
    f'{base_path}/mmpbsa_gmxmmpbsa/158/***/results/FINAL_RESULTS_MMPBSA.dat' : ((2,3,8), ('gmxmmpbsa')),
    f'{base_path}/mmpbsa_gmxmmpbsa/158/***/grinn_output/system_dry.pdb' : ((2,3,8), ('sequence')),
    }, handles={'10E8-MLA158 mod. ***' : (2,3,8)}
)

mark_resis=[pos+(240-14)-1 for pos in (34,35,36,37,63,64,65,89,90,91,92,93)]
MLA158[0].split = [(0,130),(240,350)]
MLA158[0].DSSP_assign()
MLA158[0].get_GRINN_datasets(f'{base_path}/mmpbsa_gmxmmpbsa/158/2/grinn_output')
MLA158[1].split = [(0,130),(240,350)]
MLA158[1].DSSP_assign()
MLA158[1].get_GRINN_datasets(f'{base_path}/mmpbsa_gmxmmpbsa/158/3/grinn_output')
MLA158[2].split = [(0,130),(240,350)]
MLA158[2].DSSP_assign()
MLA158[2].get_GRINN_datasets(f'{base_path}/mmpbsa_gmxmmpbsa/158/8/grinn_output')


compareMLA158 = Comparator(proteinModels=MLA158)
#model23_delta = compare158.compare_IEM(modelIndexes=[0,1])
compareMLA158.plot_IEM(modelIndexes=[0,1,2], write_best_pairs=False, mark_resis=mark_resis, dataset='total_IEM')
compareMLA158.show()


MLAWT = proteinModel.quick_load(
    {
    f'{base_path}/mmpbsa_gmxmmpbsa/WT/***/results/FINAL_RESULTS_MMPBSA.dat' : ((6,8), ('gmxmmpbsa')),
    f'{base_path}/mmpbsa_gmxmmpbsa/WT/***/grinn_output/system_dry.pdb' : ((6,8), ('sequence')),
    }, handles={'10E8-MLAWT mod. ***' : (6,8)}
)

mark_resis=[pos+(240-14)-1 for pos in (34,35,36,37,63,64,65,89,90,91,92,93)]
MLAWT[0].split = [(0,130),(240,350)]
MLAWT[0].DSSP_assign()
MLAWT[0].get_GRINN_datasets(f'{base_path}/mmpbsa_gmxmmpbsa/WT/6/grinn_output')
MLAWT[1].split = [(0,130),(240,350)]
MLAWT[1].DSSP_assign()
MLAWT[1].get_GRINN_datasets(f'{base_path}/mmpbsa_gmxmmpbsa/WT/8/grinn_output')


compareMLAWT = Comparator(proteinModels=MLAWT)
#model23_delta = compare158.compare_IEM(modelIndexes=[0,1])
compareMLAWT.plot_IEM(modelIndexes=[0,1], write_best_pairs=False, mark_resis=mark_resis, dataset='total_IEM')
compareMLAWT.show()

