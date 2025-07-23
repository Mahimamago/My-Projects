#!/usr/bin/env python
# coding: utf-8

# In[1]:


# Downloading library to create dashboard
# get_ipython().system('pip install streamlit')


# In[9]:


# Libraries
import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


# In[ ]:


# Importing Final data with model results
data = pd.read_csv('Finaldata_with_Modelresults.csv')


# In[3]:


data.head()


# In[5]:


# Adding simulated time stamp to create time summaries
base_time = pd.to_datetime('2025-01-01')
data['Simulated_Time'] = base_time + pd.to_timedelta(data.index, unit = 's')


# In[7]:


# Creating time buckets
data['Hour_Bin'] = data['Simulated_Time'].dt.floor('H')
data['12H_Bin'] = data['Simulated_Time'].dt.floor('12H')
data['24H_Bin'] = data['Simulated_Time'].dt.floor('24H')

data.head()


# #### Creating the dashboard

# In[14]:


# Starting the streamlit app
st.title("Incident Summary Dashboard")

# Sidebar timeframe selection
st.sidebar.title("TimeFrame filter")
timeframe = st.sidebar.selectbox("Select Timeframe",['Hourly', '12 Hours', '24 Hours'])

# Raw data checkbox
if st.checkbox("Show Raw data"):
    st.write(data.head())

# Selecting the appropriate bin column
if timeframe == 'Hourly':
    time_bin_col = 'Hour_Bin'
elif timeframe == '12 Hours':
    time_bin_col = '12H_Bin'
else:
    time_bin_col = '24H_Bin'


# In[20]:


# Grouping by time bin and risk label
st.subheader(f"📊 Incident Summary by Risk Label - {timeframe}")
risk_summary = data.groupby([time_bin_col, 'Risk_label']).size().unstack(fill_value = 0)
st.dataframe(risk_summary.astype(int), use_container_width=True)


# In[22]:


# Line chart of threats over time
st.subheader(f"Risk Trend over Time - {timeframe}")
fig, ax = plt.subplots(figsize=(12, 6))
risk_summary.plot(ax=ax, marker='o')
plt.xticks(rotation=45)
plt.ylabel("Number of Incidents")
plt.xlabel("Time")
plt.title(f"Trends of Risk Levels ({timeframe})")
st.pyplot(fig)


# In[26]:


# Over Risk distribution
st.subheader("Overall Risk Distribution")
fig2, ax2 = plt.subplots()
sns.countplot(data = data, x = 'Risk_label', order = data['Risk_label'].value_counts().index, palette = 'coolwarm', ax = ax2)
plt.title("Risk Label Distribution")
plt.xticks(rotation = 45)
st.pyplot(fig2)


# In[32]:


# Downloading the summary
csv = risk_summary.to_csv().encode()
st.download_button(label = "Download summary csv", data = csv, file_name = "risk_summary.csv", mime = "text/csv")


# In[30]:


# Footer
st.markdown("---")
st.markdown("Built using Streamlit")


# In[ ]:





# In[ ]:




