#!/usr/bin/env python
from netCDF4 import Dataset
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys
import os

# Usage: lcms-pda <filename> <wavelength1> <wavelength2> <wavelength3> ....
# Waters LCMS から出力したnetCDFファイルのうち最後のファイル。
# たとえばMS methodがMSscanのみ（←xxxx01.CDF）の場合は、xxxx02.CDFがAbsデータ。

arguments = sys.argv[1:]

if len(sys.argv) < 2:
    print("Usage: lcms-pda <filename> <chart_wavelength>")
    sys.exit()
else:
    file_path = sys.argv[1]
    wavelength_array = [int(wl) for wl in sys.argv[2:]]

# Load the netCDF file
nc = Dataset(file_path, mode='r')

# store netCDF data
nc_scan_acquisition_time = nc.variables['scan_acquisition_time'][:] # 秒
nc_scan_index = nc.variables['scan_index'][:] # それぞれのscanの頭の値
nc_point_count = nc.variables['point_count'][:] # 1scanあたりのデータ数
nc_mass_values = nc.variables['mass_values'][:] # m/zまたは波長の値、1scan x scan数
nc_intensity_values = nc.variables['intensity_values'][:] # データ
nc_source_file_reference = nc.getncattr('source_file_reference') #ファイル名

# Close the netCDF file
nc.close()


############################################################################################
# MS もやる！

ms_file_path = file_path[:-6]+'01.CDF'

# Load the netCDF file
msnc = Dataset(ms_file_path, mode='r')

# store netCDF data
msnc_scan_acquisition_time = msnc.variables['scan_acquisition_time'][:] # 秒
msnc_scan_index = msnc.variables['scan_index'][:] # それぞれのscanの頭の値
msnc_point_count = msnc.variables['point_count'][:] # 1scanあたりのデータ数
msnc_mass_values = msnc.variables['mass_values'][:] # m/zまたは波長の値、1scan x scan数
msnc_intensity_values = msnc.variables['intensity_values'][:] # データ
msnc_source_file_reference = msnc.getncattr('source_file_reference') #ファイル名

# Close the netCDF file
msnc.close()


############################################################################################
# Functions:
# - calculate_tic:					TIC計算
# - create_3DPDA:				3DPDAのpandas df 作製
# - export_trimmed3DPDA: 出力用にはサイズがでかすぎるのでtrim
# - get_peak_times:				クロマトグラムのピークの保持時間を得る
# - plot_chart_peak:				クロマトグラムチャートのプロット。ピーク保持時間を与えるとピーク頂点にラベル。
# - get_abs_max:					各吸収スペクトルの極大吸収の波長を得る
# - plot_spectra_subplots:	与えた保持時間の吸収スペクトルの図をプロット
# - create_msdata_df:			MS Scan の 3D df 作製
# - plot_msspectra_subplots: 与えた保持時間（吸光度保持時間+0.07 min) のマススペクトルの図をプロット

def calculate_tic(nc_scan_index,nc_intensity_values):
	"""
	Calculate total absorbance (tic) for each point (scan).
	
	Params:
		nc_scan_index
		nc_intensity_values
	
	Returns:
		np.array: TIC intensities.
	"""
	tic = []
	for i in range(len(nc_scan_index)):
	    start_index = nc_scan_index[i]
	    if i + 1 < len(nc_scan_index):
	        end_index = nc_scan_index[i + 1]
	    else:
	        end_index = len(nc_intensity_values)
	    
	    # Sum up the intensities for the current scan
	    scan_intensity_sum = np.sum(nc_intensity_values[start_index:end_index])
	    tic.append(scan_intensity_sum)
	
	# Convert TIC list to a NumPy array for easier handling
	tic = np.array(tic)
	return tic

def create_3DPDA(nc_scan_index,nc_mass_values,nc_intensity_values,time_in_minutes):
	"""
	Create 3DPDA data into dataframe.
	Params:
		nc_scan_index
		nc_mass_values
		nc_intensity_values
		time_in_minutes
	Return:
		dataframe: 3DPDA data
	"""
	spec_xaxis = []
	spec_yaxis_2D = []
	for i in range(len(nc_scan_index)):
	    start_index = nc_scan_index[i] # start_index: それぞれの　scanのfirst point
	    # end_index: それぞれのscanのlast point
	    if i + 1 < len(nc_scan_index):
	        end_index = nc_scan_index[i + 1]
	    else: # the last scan
	        end_index = len(nc_mass_values)
	    
	    # Extract the spectrum values for the current scan
	    if i == 1:
	        spec_xaxis = nc_mass_values[start_index:end_index].data
	        spec_xaxis = [int(x) for x in spec_xaxis]
	    
	    # Extract the spectrum values for the current scan
	    spec_yaxis = nc_intensity_values[start_index:end_index].data
	    spec_yaxis_2D.append(spec_yaxis)
	
	df_3D = pd.DataFrame(spec_yaxis_2D,columns=spec_xaxis,index=time_in_minutes)
	return df_3D

def export_trimmed3DPDA(df_3D):
	"""
	trim 3DPDA and save as excel sheet
	
	Parameters:
		df_3D (dataframe): large 3DPDA data
	"""
	reduced_df = df_3D.iloc[::30]
	reduced_df.to_excel(f'{outputfolder}/3DPDA.xlsx')

def get_peak_times(time_in_minutes, intensity_data, cutoff=0.05, min_time_diff=0.1):
    """
    Finds the times corresponding to peaks in a chromatogram.
    Params:
        time_in_minutes: A numpy array of time values (x-axis).
        intensity_data: A numpy array of intensity values (y-axis).
        cutoff: Minimum peak height as a percentage of the max peak.
        min_time_diff: minimum time difference.
    Returns:
        A numpy array of peak times (x-values).
    """
    max_peak_height = np.max(intensity_data)
    min_peak_height = cutoff * max_peak_height

    # Find peak indices
    peaks, _ = find_peaks(intensity_data, height=min_peak_height)

    # Get peak times
    peak_times = time_in_minutes[peaks]

    filtered_times = []
    if len(peak_times) > 0:
        filtered_times.append(peak_times[0])  # Add the first peak

        for i in range(1, len(peak_times)):
            if abs(peak_times[i] - filtered_times[-1]) >= min_time_diff:
                filtered_times.append(peak_times[i])
            elif peak_times[i] > filtered_times[-1]:
                filtered_times = filtered_times[:-1]
                filtered_times.append(peak_times[i])

    return np.array(filtered_times)

def plot_chart_peak(time_in_minutes,data,peak_times=[],title='chromatogram'):
	"""
	Plot chromatogram chart with peaks, and save as png.
	Parameters:
		time_in_minutes (np.array): x-axis.
		data (np.array): y-axis.
		peak_times (np.array): list of time of peaks to plot.
		title (string): title and file name.
	"""
	plt.figure(figsize=(10, 3))
	plt.plot(time_in_minutes, data, '-')
	if not len(peak_times) == 0: plt.plot(peak_times, data[np.isin(time_in_minutes, peak_times)], "x", color="red") #mark peaks.
	plt.title(title)
	plt.xticks(np.arange(0, 15, 1))  # Start, stop, step
	plt.xlabel('Time (min)')
	plt.grid(True, color='gray', linestyle='--', linewidth=0.5)
	plt.minorticks_on()
	plt.grid(True, which='minor', linestyle=':', linewidth='0.5', color='lightgray')
	plt.ylabel('Absorbance')
	plt.savefig(f'{outputfolder}/{title}.png')


def get_abs_max(intensity_data,x_axis,x_range=[250,600],cutoff_y=0.2):
	"""
	Finds absorbance maxima of spectrum
	Params:
		intensity_data
		x_axis (wavelength)
	"""
	selected_indices = (x_axis >= x_range[0]) & (x_axis <= x_range[1])
	selected_intensity = intensity_data[selected_indices]
	selected_x_axis = x_axis[selected_indices]
	
	max_peak_height = np.max(selected_intensity)
	min_peak_height = 1000 #cutoff_y * max_peak_height
	
	# Find abs max
	peaks, _ = find_peaks(selected_intensity, height=min_peak_height)
	# get abs max
	abs_max = selected_x_axis[peaks]

	min_abs_diff = 5
	filtered_abs = []
	if len(abs_max) > 0:
		filtered_abs.append(abs_max[0])  # Add the first peak
		
		for i in range(1, len(abs_max)):
			if abs(abs_max[i] - filtered_abs[-1]) >= min_abs_diff:
				filtered_abs.append(abs_max[i])
			elif abs_max[i] > filtered_abs[-1]:
				filtered_abs = filtered_abs[:-1]
				filtered_abs.append(abs_max[i])
	return filtered_abs

def plot_spectra_subplots(retentiontime_list, df_PDA,title='Abs Spec',cols=3):
    """
    Plots spectra as subplots with specified number of columns.

    Params:
        retentiontime_list: List of retention times to plot.
        cols: Number of columns in the subplot grid.
    """

    num_plots = len(retentiontime_list) # number of spectrum to draw
    rows = int(np.ceil(num_plots / cols))  # Calculate number of rows

    # create fig
    fig, axes = plt.subplots(rows, cols, figsize=(3 * cols, 2 * rows))  # Adjust size
    axes = axes.flatten()

    for i, time in enumerate(retentiontime_list):
        if i >= len(axes): #stop plotting if there are more subplots than data points.
            break
        
        # get spectrum
        index = np.abs(df_PDA.index.to_numpy() - time).argmin() # find closest index
        intensity = df_PDA.iloc[index].to_numpy()
        x_axis = df_PDA.columns.to_numpy()
        
        # get absmax
        absmax = get_abs_max(intensity,x_axis)
        absmax_intensity =  intensity[np.isin(x_axis, absmax)]
        print(absmax)

        ax = axes[i] # only works when there are several axes
        ax.plot(x_axis, intensity, "-")
        #if not len(absmax) == 0: ax.text(absmax, intensity[np.isin(x_axis, absmax)], "f'{absmax}'", color="red") 
        #ax.plot(absmax, intensity[np.isin(x_axis, absmax)], "x", color="red") 
        #ax.text(mz,intensity,f'{round(mz)}',fontsize=8,color='red',ha='center',va='bottom')
        ax.set_title(f'{round(time, 2)} min')
        ax.set_xlabel('Wavelength (nm)')
        ax.set_xlim(200, 700)
        ax.set_ylabel('Intensity')
        ax.set_ylim([(-1.2)*np.abs(np.min(intensity)), np.max(intensity[100:])*1.2])
        #ax.set_yticks([0,20000,40000])
        #ax.set_yticklabels([0,2,4])
        for wl, value in zip(absmax,absmax_intensity):
        	ax.text(wl,value,f'{round(wl)}',fontsize=8,color='red',ha='center')


    # Remove extra subplots if needed
    if num_plots < len(axes):
        for j in range(num_plots, len(axes)):
            fig.delaxes(axes[j])

    plt.tight_layout()
    plt.savefig(f'{outputfolder}/{title}.png')


############################################################################################
# functions for MS


def create_msdata_df(nc_scan_index,nc_mass_values,nc_intensity_values,msnc_scan_acquisition_time):
	"""
	Create DataFrame of mass spec data
	Params:
		nc_scan_index
		nc_mass_values
		nc_intensity_values
		nc_scan_acquisition_time
	Return:
		dataframe: 3DPDA data
	"""
	time_in_minutes = msnc_scan_acquisition_time / 60
	mz_list = []
	intensity_list = []
	for i in range(len(nc_scan_index)):
	    start_index = nc_scan_index[i] # start_index: それぞれの　scanのfirst point
	    # end_index: それぞれのscanのlast point
	    if i + 1 < len(nc_scan_index):
	        end_index = nc_scan_index[i + 1]
	    else: # the last scan
	        end_index = len(nc_mass_values)
	    
	    # Extract the spectrum values for the current scan
	    mz = nc_mass_values[start_index:end_index].data
	    mz_list.append(mz)
	    
	    # Extract the spectrum values for the current scan
	    intensity = nc_intensity_values[start_index:end_index].data
	    intensity_list.append(intensity)
	
	df_msdata = pd.DataFrame({'m/z': mz_list, 'Intensity': intensity_list},index=time_in_minutes)
	return df_msdata

def plot_msspectra_subplots(retentiontime_list, df_msdata,title='MS Spec',percent_cutoff='50',cols=3):
    """
    Plots spectra as subplots with specified number of columns.

    Params:
        retentiontime_list: List of retention times to plot.
        cols: Number of columns in the subplot grid.
    """

    num_plots = len(retentiontime_list) # number of spectrum to draw
    rows = int(np.ceil(num_plots / cols))  # Calculate number of rows

    # create fig
    fig, axes = plt.subplots(rows, cols, figsize=(3 * cols, 2 * rows))  # Adjust size
    axes = axes.flatten()

    for i, time in enumerate(retentiontime_list):
        if i >= len(axes): #stop plotting if there are more subplots than data points.
            break
        
        mstime = time + 0.07
        
        # get spectrum
        index = np.abs(df_msdata.index.to_numpy() - mstime).argmin() # find closest index
        msspec = df_msdata.iloc[index]
        mz_list = msspec['m/z']
        intensity_list = msspec['Intensity']

        ax = axes[i] # only works when there are several axes
        ax.bar(mz_list,intensity_list,width=5,align='center')
        ax.set_title(f'{round(time, 2)}+0.07 min')
        ax.set_xlabel('m/z')
        #ax.set_xlim(200, 700)
        ax.set_ylabel('Intensity')
        ax.set_ylim([(-1.2)*np.abs(np.min(intensity_list)), np.max(intensity_list)*1.3])
        #ax.set_yticks([0,20000,40000])
        #ax.set_yticklabels([0,2,4])

        # Define the threshold as a percentage of the maximum intensity
        percentage_threshold = int(percent_cutoff)  # For example, 50%
        max_intensity = max(intensity_list)
        threshold = (percentage_threshold / 100) * max_intensity
        for mz, intensity in zip(mz_list,intensity_list):
            if intensity > threshold:
                ax.text(mz,intensity,f'{round(mz)}',fontsize=8,color='red',ha='center',va='bottom')

    # Remove extra subplots if needed
    if num_plots < len(axes):
        for j in range(num_plots, len(axes)):
            fig.delaxes(axes[j])

    plt.tight_layout()
    plt.savefig(f'{outputfolder}/{title}.png')



############################################################################################

if __name__ == "__main__":
	
	# Convert time from seconds to minutes
	time_in_minutes = nc_scan_acquisition_time / 60
	
	#create output folder
	outputfolder = nc_source_file_reference+'-output'
	if not os.path.exists(outputfolder):
	    os.makedirs(outputfolder)

	############### create PDA and export as excel ###############
	df_PDA = create_3DPDA(nc_scan_index,nc_mass_values,nc_intensity_values,time_in_minutes)
	# export_trimmed3DPDA(df_PDA)
	
	################## TIC ##################
	### calculate tic 
	#tic_data = calculate_tic(nc_scan_index,nc_intensity_values)
	#plot_chart_peak(time_in_minutes,tic_data,title='TIC')
	
	### get retention times of peaks
	#peak_times_tic = get_peak_times(time_in_minutes,tic_data,cutoff=0.05)
	
	### plot TIC chart with peaks
	#plot_chart_peak(time_in_minutes,tic_data,peak_times_tic,f"TIC - Chart")
	
	### plot spectra
	#plot_spectra_subplots(peak_times_tic,df_PDA,f"TIC - AbsSpec")

	################## MS ##################
	df_msdata = create_msdata_df(msnc_scan_index,msnc_mass_values,msnc_intensity_values,msnc_scan_acquisition_time)
	
	############### selected wavelength ###############
	for wavelength in wavelength_array:
		# select wavelength (nm) value you're interested in, and get the data
		sw_data = df_PDA[wavelength].to_numpy()
		
		# get retention times of peaks
		peak_times = get_peak_times(time_in_minutes,sw_data,cutoff=0.1,min_time_diff=0.05)
		
		# plot chromatogram chart with peaks
		plot_chart_peak(time_in_minutes,sw_data,peak_times,f"{wavelength} nm - Chart")
		
		# plot spectra
		plot_spectra_subplots(peak_times,df_PDA,f"{wavelength} nm - AbsSpec")
	
		# plot MS spectra
		plot_msspectra_subplots(peak_times, df_msdata,title=f'{wavelength} - MSSpec',percent_cutoff='40',cols=3)
	
