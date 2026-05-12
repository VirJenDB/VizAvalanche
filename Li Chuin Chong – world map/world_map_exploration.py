from pathlib import Path
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
import numpy as np

CSV_PATH = Path("data/collected_data.csv")

df = pd.read_csv(CSV_PATH)
print(f"Loaded: {CSV_PATH}")
print(f"Rows: {df.shape[0]:,}")
print(f"Columns: {df.shape[1]:,}")

df.head()

# inspect the column
df["submitter_country"].isna().sum() # 192
df["submitter_country"].value_counts()

country_counts = (
    df["submitter_country"]
    .value_counts()
    .reset_index()
)

country_counts.columns = ["country", "count"]
country_counts

# log count
import numpy as np

country_counts["log_count"] = np.log10(country_counts["count"] + 1)
country_counts.head()

max_count = country_counts["count"].max()

tick_counts = [1, 10, 100, 1_000, 10_000, 100_000, 1_000_000,]
if max_count >= 5_000_000:
    tick_counts.append(5_000_000)
if max_count >= 10_000_000:
    tick_counts.append(10_000_000)
tick_counts = [x for x in tick_counts if x <= max_count or x == 1]

tick_vals = [np.log10(x + 1) for x in tick_counts]

def format_count_label(x):
    if x >= 1_000_000:
        return f"{x // 1_000_000}M"
    elif x >= 1_000:
        return f"{x // 1_000}k"
    else:
        return str(x)
tick_text = [format_count_label(x) for x in tick_counts]

fig_map = go.Figure(
    data=go.Choropleth(
        locations=country_counts["country"],
        z=country_counts["log_count"],
        locationmode="country names",
        text=country_counts["country"],
        customdata=country_counts[["count"]],
        colorbar=dict(
            title="Number of records",
            tickvals=tick_vals,
            ticktext=tick_text,
        ),
        hovertemplate=(
            "<b>%{text}</b><br>"
            "Records: %{customdata[0]:,}<br>"
            "<extra></extra>"
        ),
        marker_line_color="white",
        marker_line_width=0.5,
    )
)

fig_map.update_layout(
    title="World map of submitter countries (log-scaled by record count)",
    geo=dict(
        showframe=False,
        showcoastlines=True,
        projection_type="natural earth",
    ),
    margin=dict(l=0, r=0, t=60, b=0),
)

fig_map.show()
