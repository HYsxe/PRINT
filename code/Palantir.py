import h5py
import numpy as np
import palantir
import pandas as pd

# Transfer data from R to python
h5_path = pd.read_csv("palantir_args.txt", header = None).iloc[0, 0]
print(h5_path)
hf = h5py.File(h5_path, 'r')
cell_cisTopic = np.transpose(hf.get('cellCisTopic'))
barcodes = hf.get('barcodes')
root_barcode = hf.get('rootCellBarcode')
barcodes = [bc.decode('ascii') for bc in barcodes]
root_barcode = [bc.decode('ascii') for bc in root_barcode][0]
hf.close()

# Create data frame of cisTopic embedding
embedding_df = pd.DataFrame(data = cell_cisTopic,
                            index = barcodes)

# Run diffusion maps
dm_res = palantir.utils.run_diffusion_maps(embedding_df, n_components=5)
ms_data = palantir.utils.determine_multiscale_space(dm_res)

# Calculate pseudotime from designated root cell
pr_res = palantir.core.run_palantir(ms_data, 
                                    root_barcode, 
                                    num_waypoints=1000,
                                    use_early_cell_as_start=True)

# Save results
hf = h5py.File(h5_path, 'w')
hf.create_dataset("pseudotime", data = pr_res.pseudotime.values)
hf.create_dataset("lineageProbs", data = pr_res.branch_probs)
hf.close()
