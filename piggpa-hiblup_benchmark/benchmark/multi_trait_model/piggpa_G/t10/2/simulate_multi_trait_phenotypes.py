import numpy as np
import pandas as pd
import os

np.random.seed(42)

output_dir = '/public/share/likui/hanyu/testdata/In-silico-data/t10/2'
os.makedirs(output_dir, exist_ok=True)

original_pheno = pd.read_csv('/public/share/likui/hanyu/testdata/In-silico-data/phenotypes.txt', sep='\t')

n_individuals = len(original_pheno)

ids = original_pheno['ID'].values
base_phenotype = original_pheno['Phenotype'].values
base_breeding_value = original_pheno['BreedingValue'].values

sex = np.random.choice(['M', 'F'], size=n_individuals)
season = np.random.choice(['Spring', 'Summer', 'Autumn', 'Winter'], size=n_individuals)

sex_effect_t1 = np.where(sex == 'M', 5.0, 0.0)
sex_effect_t2 = np.where(sex == 'M', 8.0, 0.0)
sex_effect_t3 = np.where(sex == 'M', -1.0, 0.0)

season_effects_t1 = {'Spring': 100, 'Summer': 105, 'Autumn': 95, 'Winter': 90}
season_effects_t2 = {'Spring': 110, 'Summer': 120, 'Autumn': 100, 'Winter': 95}
season_effects_t3 = {'Spring': 15, 'Summer': 18, 'Autumn': 14, 'Winter': 12}

season_effect_t1 = np.array([season_effects_t1[s] for s in season])
season_effect_t2 = np.array([season_effects_t2[s] for s in season])
season_effect_t3 = np.array([season_effects_t3[s] for s in season])

base_t1 = 100 + base_phenotype * 10 + sex_effect_t1 + season_effect_t1
base_t2 = 110 + base_phenotype * 8 + base_breeding_value * 5 + sex_effect_t2 + season_effect_t2
base_t3 = 15 + base_phenotype * 2 + sex_effect_t3 + season_effect_t3

noise_sd_t1 = 3.0
noise_sd_t2 = 5.0
noise_sd_t3 = 1.5

T1 = base_t1 + np.random.normal(0, noise_sd_t1, n_individuals)
T2 = base_t2 + np.random.normal(0, noise_sd_t2, n_individuals)
T3 = base_t3 + np.random.normal(0, noise_sd_t3, n_individuals)

T1 = np.round(T1, 1)
T2 = np.round(T2, 1)
T3 = np.round(T3, 1)

result = pd.DataFrame({
    'ID': ids,
    'T1': T1,
    'T2': T2,
    'T3': T3,
    'sex': sex,
    'season': season
})

output_file = os.path.join(output_dir, 'simulated_phenotypes_multi_trait.txt')
result.to_csv(output_file, sep='\t', index=False)

print(f"多性状模拟表型文件已保存到: {output_file}")
print(f"总个体数: {n_individuals}")
print(f"\n前10行数据预览:")
print(result.head(10).to_string(index=False))
print(f"\n数据统计:")
print(f"性别分布: M={sum(sex=='M')}, F={sum(sex=='F')}")
print(f"季节分布: Spring={sum(season=='Spring')}, Summer={sum(season=='Summer')}, Autumn={sum(season=='Autumn')}, Winter={sum(season=='Winter')}")
print(f"\n性状统计:")
print(f"T1: mean={T1.mean():.2f}, std={T1.std():.2f}, min={T1.min():.2f}, max={T1.max():.2f}")
print(f"T2: mean={T2.mean():.2f}, std={T2.std():.2f}, min={T2.min():.2f}, max={T2.max():.2f}")
print(f"T3: mean={T3.mean():.2f}, std={T3.std():.2f}, min={T3.min():.2f}, max={T3.max():.2f}")

corr_t1_t2 = np.corrcoef(T1, T2)[0, 1]
corr_t1_t3 = np.corrcoef(T1, T3)[0, 1]
corr_t2_t3 = np.corrcoef(T2, T3)[0, 1]
print(f"\n性状间相关性:")
print(f"T1-T2: r={corr_t1_t2:.3f}")
print(f"T1-T3: r={corr_t1_t3:.3f}")
print(f"T2-T3: r={corr_t2_t3:.3f}")
