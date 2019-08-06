from grit.lib.multiprocessing_utils import run_in_parallel
all_args = []
for ct in process_data.iter_ct('all'):
    all_args.append([ct, ['mean', 'max', 'min', 'q25']])
run_in_parallel(8, , all_args)
