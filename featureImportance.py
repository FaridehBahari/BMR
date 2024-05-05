import numpy as np
import pandas as pd

# Example data dimensions
num_samples = 1000
num_genomic_epigenomic_features = 50
num_llm_features = 30

# Generate example genomic/epigenomic features
X_genomic_epigenomic = np.random.randn(num_samples, num_genomic_epigenomic_features)

# Generate example LLM features
X_llm = np.random.randn(num_samples, num_llm_features)

# Generate example target variable
y = np.random.rand(num_samples)

# Create DataFrames for genomic/epigenomic and LLM features
genomic_epigenomic_columns = ['GenomicEpigenomicFeature{}'.format(i) for i in range(1, num_genomic_epigenomic_features + 1)]
llm_columns = ['LLMFeature{}'.format(i) for i in range(1, num_llm_features + 1)]

df_genomic_epigenomic = pd.DataFrame(X_genomic_epigenomic, columns=genomic_epigenomic_columns)
df_llm = pd.DataFrame(X_llm, columns=llm_columns)

# Concatenate features and target variable into a single DataFrame
df_genomic_epigenomic['Target'] = y
df_llm['Target'] = y

# Display sample of the data
print("Genomic/Epigenomic Features:")
print(df_genomic_epigenomic.head())

print("\nLLM Features:")
print(df_llm.head())


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split

# Assuming you have separate datasets for genomic/epigenomic features and LLM features
# X_genomic_epigenomic and X_llm represent feature matrices
# y represents the target variable

# Split data into training and testing sets for genomic/epigenomic features
X_train_genomic_epigenomic, X_test_genomic_epigenomic, y_train, y_test = train_test_split(
    X_genomic_epigenomic, y, test_size=0.2, random_state=42)

# Train a random forest regressor model on genomic/epigenomic features
rf_model_genomic_epigenomic = RandomForestRegressor(n_estimators=100, random_state=42)
rf_model_genomic_epigenomic.fit(X_train_genomic_epigenomic, y_train)

# Calculate feature importances for genomic/epigenomic model
feature_importances_genomic_epigenomic = rf_model_genomic_epigenomic.feature_importances_

# Split data into training and testing sets for LLM features
X_train_llm, X_test_llm, y_train, y_test = train_test_split(
    X_llm, y, test_size=0.2, random_state=42)

# Train a random forest regressor model on LLM features
rf_model_llm = RandomForestRegressor(n_estimators=100, random_state=42)
rf_model_llm.fit(X_train_llm, y_train)

# Calculate feature importances for LLM model
feature_importances_llm = rf_model_llm.feature_importances_

# Plot feature importances for both models
plt.figure(figsize=(10, 6))
plt.bar(range(len(feature_importances_genomic_epigenomic)), feature_importances_genomic_epigenomic,
        color='skyblue', label='Genomic/Epigenomic Features')
plt.bar(range(len(feature_importances_llm)), feature_importances_llm,
        color='orange', alpha=0.7, label='LLM Features')
plt.title('Feature Importance Comparison')
plt.xlabel('Feature Index')
plt.ylabel('Feature Importance')
plt.legend()
plt.show()


# Sort feature importances arrays in decreasing order
sorted_indices_genomic_epigenomic = np.argsort(feature_importances_genomic_epigenomic)[::-1]
sorted_indices_llm = np.argsort(feature_importances_llm)[::-1]

sorted_feature_importances_genomic_epigenomic = feature_importances_genomic_epigenomic[sorted_indices_genomic_epigenomic]
sorted_feature_importances_llm = feature_importances_llm[sorted_indices_llm]

# Plot sorted feature importances for both models
plt.figure(figsize=(10, 6))
plt.bar(range(len(sorted_feature_importances_genomic_epigenomic)), sorted_feature_importances_genomic_epigenomic,
        color='skyblue', label='Genomic/Epigenomic Features')
plt.bar(range(len(sorted_feature_importances_llm)), sorted_feature_importances_llm,
        color='orange', alpha=0.7, label='LLM Features')
plt.title('Feature Importance Comparison')
plt.xlabel('Feature Index')
plt.ylabel('Feature Importance')
plt.legend()
plt.show()

##################################################
##################################################
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split

# Example data dimensions
num_samples = 1000
num_genomic_epigenomic_features = 50
num_llm_features = 30

# Generate example genomic/epigenomic features
X_genomic_epigenomic = np.random.randn(num_samples, num_genomic_epigenomic_features)

# Generate example LLM features
X_llm = np.random.randn(num_samples, num_llm_features)

# Generate example target variable
y = np.random.rand(num_samples)

# Create DataFrames for genomic/epigenomic and LLM features
genomic_epigenomic_columns = ['GenomicEpigenomicFeature{}'.format(i) for i in range(1, num_genomic_epigenomic_features + 1)]
llm_columns = ['LLMFeature{}'.format(i) for i in range(1, num_llm_features + 1)]

df_genomic_epigenomic = pd.DataFrame(X_genomic_epigenomic, columns=genomic_epigenomic_columns)
df_llm = pd.DataFrame(X_llm, columns=llm_columns)

# Concatenate features and target variable into a single DataFrame
df_genomic_epigenomic['Target'] = y
df_llm['Target'] = y

# Load feature groups from Excel file
feature_groups_df = pd.read_excel('../external/BMR/rawInput/41467_2019_13929_MOESM5_ESM.xlsx')  # Update with your file path
feature_groups = feature_groups_df.groupby('Group Name')['Feature Name'].apply(list).to_dict()

# Calculate feature importances for genomic/epigenomic features
rf_model_genomic_epigenomic = RandomForestRegressor(n_estimators=100, random_state=42)
rf_model_genomic_epigenomic.fit(X_genomic_epigenomic, y)

feature_importances_genomic_epigenomic = rf_model_genomic_epigenomic.feature_importances_

# Calculate feature importances for LLM features
rf_model_llm = RandomForestRegressor(n_estimators=100, random_state=42)
rf_model_llm.fit(X_llm, y)

feature_importances_llm = rf_model_llm.feature_importances_

# Aggregate feature importance for genomic/epigenomic features by group
feature_importance_grouped = {}
for group, features in feature_groups.items():
    indices = [genomic_epigenomic_columns.index(feature) for feature in features if feature in genomic_epigenomic_columns]
    group_importance = np.sum(feature_importances_genomic_epigenomic[indices])
    feature_importance_grouped[group] = group_importance

# Plot feature importances for both models
plt.figure(figsize=(10, 6))
plt.bar(range(len(feature_importance_grouped)), feature_importance_grouped.values(),
        color='skyblue', label='Genomic/Epigenomic Features Grouped')
plt.bar(range(len(feature_importances_llm)), feature_importances_llm,
        color='orange', alpha=0.7, label='LLM Features')
plt.xticks(range(len(feature_importance_grouped)), feature_importance_grouped.keys(), rotation=45, ha='right')
plt.title('Feature Importance Comparison')
plt.xlabel('Feature Group')
plt.ylabel('Feature Importance')
plt.legend()
plt.tight_layout()
plt.show()
