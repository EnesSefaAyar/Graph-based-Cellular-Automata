from joblib import Parallel, delayed
import main_barabasi as m

res = Parallel(n_jobs=-1)(delayed(m.main)(i, (i+1)*200) for i in range(0, 100))

