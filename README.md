# waters-lcms-netCDF-analysis
* analyze Waters LCMS netCDF files to visualize LC charts and Abs/MS spectra
* Waters LC-PDA-MS の netCDFファイルを使って、目的吸収波長のクロマトグラムを描き、ピークの保持時間の吸収スペクトルとそれに対応するマススペクトルを描画

## Usage
`python this_script.py filename.CDF wavelength1 wavelength2 wavelength3 ....`
* filename.CDF: PDAデータのnetCDF ファイル = Waters LCMS から出力したnetCDFファイルのうち最後のもの。
* たとえばMS methodがMSscanのみ（←xxxx01.CDF）の場合は、xxxx02.CDFがAbsデータ。
