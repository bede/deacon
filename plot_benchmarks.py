#!/usr/bin/env python3
"""
Plot benchmark results from benchmark_results.tsv
Shows max_bp_per_second for lr and sr separately in chronological order
"""

import pandas as pd
import altair as alt
import subprocess

# Read the TSV file
df = pd.read_csv('benchmark_results.tsv', sep='\t')

# Convert bp/s to Mbp/s
df['max_Mbps'] = df['max_bp_per_second'] / 1_000_000

# Get commit dates from git history
commits = df['commit'].unique().tolist()
commit_dates = {}
for commit in commits:
    git_cmd = ['git', 'show', '-s', '--format=%ci', commit]
    git_output = subprocess.run(git_cmd, capture_output=True, text=True)
    date_str = git_output.stdout.strip().split()[0]  # Get just YYYY-MM-DD
    commit_dates[commit] = date_str

# Add dates to dataframe
df['date'] = df['commit'].map(commit_dates)
df['date'] = pd.to_datetime(df['date'])

# Sort by date
df = df.sort_values('date')

# Create separate dataframes for lr and sr
df_lr = df[df['type'] == 'lr'].copy()
df_sr = df[df['type'] == 'sr'].copy()

# Create the plot with lines and points
base = alt.Chart(df).encode(
    x=alt.X('date:T',
            title='Date',
            axis=alt.Axis(format='%Y-%m-%d')),
    y=alt.Y('max_Mbps:Q',
            title='Throughput (Mbp/s)',
            scale=alt.Scale(zero=False)),
    color=alt.Color('type:N',
                    title='Test',
                    scale=alt.Scale(domain=['lr', 'sr'],
                                   range=['#1f77b4', '#ff7f0e']),
                    legend=alt.Legend(
                        orient='top-right',
                        labelExpr="datum.label == 'lr' ? 'Single reads' : 'Paired reads'"
                    ))
)

lines = base.mark_line()
points = base.mark_circle(size=100).encode(
    tooltip=[
        alt.Tooltip('date:T', title='Date', format='%Y-%m-%d'),
        alt.Tooltip('commit:N', title='Commit'),
        alt.Tooltip('type:N', title='Type'),
        alt.Tooltip('max_Mbps:Q', title='Throughput (Mbp/s)', format='.2f'),
        alt.Tooltip('trial1_bp_per_second:Q', title='Trial 1 (bp/s)', format=','),
        alt.Tooltip('trial2_bp_per_second:Q', title='Trial 2 (bp/s)', format=',')
    ]
)

chart = (lines + points).properties(
    width=800,
    height=500,
    title='Deacon M1 Pro throughput by commit (rsviruses17900; single: fasta, paired: fq.zst)'
).interactive()

# Save to HTML
chart.save('benchmark_plot.html')
print("Plot saved to benchmark_plot.html")

# Also display statistics
print("\nSummary Statistics:")
print("\nSingle Reads (lr):")
print(f"  Mean: {df_lr['max_Mbps'].mean():.2f} Mbp/s")
print(f"  Min:  {df_lr['max_Mbps'].min():.2f} Mbp/s (commit {df_lr.loc[df_lr['max_Mbps'].idxmin(), 'commit']})")
print(f"  Max:  {df_lr['max_Mbps'].max():.2f} Mbp/s (commit {df_lr.loc[df_lr['max_Mbps'].idxmax(), 'commit']})")

print("\nPaired Reads (sr):")
print(f"  Mean: {df_sr['max_Mbps'].mean():.2f} Mbp/s")
print(f"  Min:  {df_sr['max_Mbps'].min():.2f} Mbp/s (commit {df_sr.loc[df_sr['max_Mbps'].idxmin(), 'commit']})")
print(f"  Max:  {df_sr['max_Mbps'].max():.2f} Mbp/s (commit {df_sr.loc[df_sr['max_Mbps'].idxmax(), 'commit']})")
