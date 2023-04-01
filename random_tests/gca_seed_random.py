from joblib import Parallel, delayed
import main_seed_random as m

cell_line = "ACH-000001"

res = Parallel(n_jobs=-1)(delayed(m.main)(cell_line, i) for i in range(100))

