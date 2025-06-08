import json
import numpy as np
import pandas as pd

# Remove invalid subcarriers
CSI_VAID_SUBCARRIER_INTERVAL = 1

csi_vaid_subcarrier_index = []
# LLTF: 52
csi_vaid_subcarrier_index += [i for i in range(6, 32, CSI_VAID_SUBCARRIER_INTERVAL)]     # 26  red
csi_vaid_subcarrier_index += [i for i in range(33, 59, CSI_VAID_SUBCARRIER_INTERVAL)]    # 26  green
CSI_DATA_LLFT_COLUMNS = len(csi_vaid_subcarrier_index)
# HT-LFT: 56 + 56
csi_vaid_subcarrier_index += [i for i in range(66, 94, CSI_VAID_SUBCARRIER_INTERVAL)]    # 28  blue
csi_vaid_subcarrier_index += [i for i in range(95, 123, CSI_VAID_SUBCARRIER_INTERVAL)]   # 28  White

CSI_DATA_INDEX = 200  # buffer size
CSI_DATA_COLUMNS = len(csi_vaid_subcarrier_index)

def raw2csi(csi_raw_data):
    print(csi_raw_data.shape)
    CSI_DATA_INDEX = len(csi_raw_data)

    csi_data_array = np.zeros(
        [CSI_DATA_INDEX, CSI_DATA_COLUMNS], dtype=np.complex64)
    
    # 将CSI_DATA整数两两组合恢复成复数形式
    for row in range(CSI_DATA_INDEX):
        for i in range(CSI_DATA_COLUMNS):
            csi_data_array[row][i] = complex(csi_raw_data[row,csi_vaid_subcarrier_index[i] * 2+1],
                                            csi_raw_data[row,csi_vaid_subcarrier_index[i] * 2])
            
    return csi_data_array


def read_csv(file_path):
    df = pd.read_csv(file_path)
    raw_data = df['data'].values
    raw_timestamp = df['local_timestamp'].values

    timestamp = []
    for i in range(len(raw_timestamp)):
        timestamp.append((raw_timestamp[i]-raw_timestamp[0])/1e6)
    timestamp = np.array(timestamp)
    # print(timestamp[:10])

    csi = []
    for i in range(len(raw_data)):
        csi.append(json.loads(raw_data[i]))

    data = raw2csi(np.array(csi))
    return data,timestamp

if __name__ == '__main__':
    data,time = read_csv("data/room_kfw/go/2023-11-28_20-05-47-773_104_1.csv")
    print(data.shape)
