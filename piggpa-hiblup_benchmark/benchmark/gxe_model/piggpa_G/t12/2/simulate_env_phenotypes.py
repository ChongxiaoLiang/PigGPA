import numpy as np
import pandas as pd
import os

np.random.seed(42)

output_dir = '/public/share/likui/hanyu/testdata/In-silico-data/t12/2'
os.makedirs(output_dir, exist_ok=True)

original_pheno = pd.read_csv('/public/share/likui/hanyu/testdata/In-silico-data/phenotypes.txt', sep='\t')

n_individuals = len(original_pheno)

ids = original_pheno['ID'].values
base_phenotype = original_pheno['Phenotype'].values
base_breeding_value = original_pheno['BreedingValue'].values

temp = np.random.uniform(20, 35, n_individuals)
humidity = np.random.uniform(50, 90, n_individuals)

temp_effect = (temp - 27.5) * 2.0
humidity_effect = (humidity - 70) * 0.5

base_trait = 100 + base_phenotype * 10 + base_breeding_value * 5

trait = base_trait + temp_effect + humidity_effect + np.random.normal(0, 5, n_individuals)

temp = np.round(temp, 1)
humidity = np.round(humidity, 1)
trait = np.round(trait, 1)

result = pd.DataFrame({
    'ID': ids,
    'Temp': temp,
    'Humidity': humidity,
    'Trait': trait
})

output_file = os.path.join(output_dir, 'simulated_phenotypes_env.txt')
result.to_csv(output_file, sep='\t', index=False)

print(f"模拟表型文件已保存到: {output_file}")
print(f"总个体数: {n_individuals}")
print(f"\n前10行数据预览:")
print(result.head(10).to_string(index=False))
print(f"\n环境变量统计:")
print(f"温度: mean={temp.mean():.2f}, std={temp.std():.2f}, min={temp.min():.2f}, max={temp.max():.2f}")
print(f"湿度: mean={humidity.mean():.2f}, std={humidity.std():.2f}, min={humidity.min():.2f}, max={humidity.max():.2f}")
print(f"\n性状统计:")
print(f"Trait: mean={trait.mean():.2f}, std={trait.std():.2f}, min={trait.min():.2f}, max={trait.max():.2f}")
print(f"\n相关性:")
print(f"Temp-Trait: r={np.corrcoef(temp, trait)[0,1]:.3f}")
print(f"Humidity-Trait: r={np.corrcoef(humidity, trait)[0,1]:.3f}")
