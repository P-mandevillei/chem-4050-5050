from compute_work_isothermal import compute_work_isothermal
from compute_work_adiabatic import compute_work_adiabatic
from matplotlib import pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np

n = 1 # mol
temp_i = 300 # K
v_i = 0.1 # m^3
gamma = 1.4

vols = np.linspace(v_i, 3*v_i, 100)
works = pd.DataFrame({'final_volume': vols})
works['isothermal'] = works['final_volume'].apply(lambda v_f: compute_work_isothermal(v_i, v_f, n, temp_i))
works['adiabatic'] = works['final_volume'].apply(lambda v_f: compute_work_adiabatic(v_i, v_f, gamma, n, temp_i))
works.rename(columns={'final_volume': 'final_volume_m^3', 'isothermal': 'isothermal_J', 'adiabatic': 'adiabatic_J'}).to_csv('work.csv')
print("Saved computed works for both isothermal and adiabatic processes to work.csv")

plt_works = works.melt(
    id_vars='final_volume',
    value_vars=['isothermal', 'adiabatic'],
    var_name='condition',
    value_name='work'
)
fig, ax = plt.subplots()
fig.set_size_inches(6,4)
sns.lineplot(plt_works, x='final_volume', y='work', hue='condition')
ax.set_ylabel('Work (J)')
ax.set_xlabel('Final Volume ($m^3$)')
ax.set_title('Work done by an ideal gas in isothermal and adiabatic expansion')
fig.savefig('work.png', dpi=210)
print("Saved plot to work.png")

with open('work_discussion.md', 'w') as file:
	file.write(r"""The work done by an ideal gas on the surroundings in an isothermal expansion is always greater than the work done in an adiabatic process for the same volume change. This makes chemical sense. $\Delta U = w + q$, in the adiabatic process, $q=0$, $w<0$, so $\Delta U<0$. Since $\Delta U = n C_v \Delta T$ for ideal gas, $\Delta T<0$. The ideal gas law states that $PV=nRT$, so if $T$ decreases, $P$ must also decrease. Therefore for the adiabatic process, the pressure integrand is also less than that in the isothermal process.""")
print("Saved discussion to work_discussion.md")