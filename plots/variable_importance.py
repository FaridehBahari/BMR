import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import pickle
import xgboost as xgb

# Load the GBM model
path = '../external/BMR/output/with_RepliSeq_HiC/bin_size_effect/var_size/GBM/GBM_model.pkl'
with open(path, 'rb') as file:
    model = pickle.load(file)

# After training the model (assuming your model is named 'model')
importance = model.get_score(importance_type='weight')  # Can also use 'gain' or 'cover'

# Sort the importance dictionary by value (importance score) in descending order
sorted_importance = {k: v for k, v in sorted(importance.items(), key=lambda item: item[1], reverse=True)}

# Print the top 50 features
top_50 = dict(list(sorted_importance.items())[:50])
for feature, score in top_50.items():
    print(f"{feature}: {score}")

# Load feature groups
with open('../external/procInput/ftrs_dict.pickle', 'rb') as file:
    feature_groups = pickle.load(file)

# Filter out specific keys
exclude_keys = ['DNA_methylation', 'APOBEC']
feature_groups_filtered = {key: value for key, value in feature_groups.items() if key not in exclude_keys}

# Define the colors
group_colors = {
    'Epigenetic_mark': "#66C2A5",
    'RNA_expression': "#3288BD",
    'Replication_timing': "#43589F",
    'HiC': "#D53E4F",
    'DNA_accessibility': "#ABDDA4",
    'nucleotide content': "#D4B9DA",
    'conservation': "#D1E5F0"
}

# Map features to their group colors
feature_color_map = {}
for group, features in feature_groups_filtered.items():
    for feature in features:
        feature_color_map[feature] = group_colors[group]

# Set the font to DejaVu Serif
plt.rcParams['font.family'] = 'DejaVu Serif'

# Set the figure size
plt.figure(figsize=(10, 10))

# Apply colors by mapping each feature to its group color
colors = [feature_color_map.get(feature, 'grey') for feature in top_50.keys()]  # default to 'grey' if not found
bars = plt.barh(range(len(top_50)), list(top_50.values()), align='center', color=colors)
plt.yticks(range(len(top_50)), list(top_50.keys()), fontsize=13)  # Increase font size for y-ticks
plt.xlabel('Importance Score', fontsize=16)  # Increase font size for x-label

# plt.title('Top 50 Feature Importance', fontsize=18)  # Increase font size for title
plt.gca().invert_yaxis()  # Invert y-axis to have the highest importance at the top

# Format legend labels
def format_legend_label(label):
    label = label.replace('_', ' ').title()
    label = label.replace('Dna', 'DNA').replace('Rna', 'RNA')
    return label

# Create a legend for the colors
legend_elements = [Patch(facecolor=color, edgecolor='w', label=format_legend_label(group)) for group, color in group_colors.items()]
plt.legend(handles=legend_elements, loc='lower right', fontsize=16)  # Increase font size for legend

# Adjust layout to ensure all labels fit
plt.tight_layout()

# Save the plot
plt.savefig('../external/BMR/plots/top_50_variable_importance.png')
plt.show()
