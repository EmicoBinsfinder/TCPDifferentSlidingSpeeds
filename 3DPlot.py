xdata = [300, 400, 500, 600, 700]
ydata = [1, 2, 3, 4, 5]

Temperatures = ["300K", "400K", "500K", "600K", "700K"]
df = pd.DataFrame(columns=[Temperatures])
df.loc[len(df)] = Average_Shear_Stress_List_1GPa
df.loc[len(df)] = Average_Shear_Stress_List_2GPa
df.loc[len(df)] = Average_Shear_Stress_List_3GPa
df.loc[len(df)] = Average_Shear_Stress_List_4GPa
df.loc[len(df)] = Average_Shear_Stress_List_5GPa

#print(df)

df = df.T
#print(df)

Average_Shear_Stress_List_1GPa = sorted(df.iloc[0])
Average_Shear_Stress_List_2GPa = sorted(df.iloc[1])
Average_Shear_Stress_List_3GPa = sorted(df.iloc[2])
Average_Shear_Stress_List_4GPa = sorted(df.iloc[3])
Average_Shear_Stress_List_5GPa = sorted(df.iloc[4])

#print(Average_Shear_Stress_List_1GPa)

def function(data, a, b, c):
    x = data[0]
    y = data[1]
    return a * (x**b) * (y**c) #TODO change fitting function

x_data = []
y_data = []
z_data = []

data = [[300, 0.328462898923518,  0.6061013664596271], [300, 0.3351728917540594, 3.011891925465839], [300, 0.3748952949778793, 4.012132919254659], [300, 0.3821062799760507, 7.509253416149069], [300, 0.4078641680450279, 9.997445962732918],
        [400, 0.363257762286299, 1.853009999999999], [400, 0.3820543502176791, 4.6494], [400, 0.42288773849385886, 9.3663], [400, 0.42288773849385886, 13.273200000000001], [400, 0.4291421178005783, 18.4923],
        [500, 0.4376893408917189, 3.2911218274111675], [500, 0.4509062860550436, 9.685736040609136], [500, 0.4662144565589674, 12.320154822335025], [500, 0.47616665712067596, 21.711769035532996], [500, 0.47756143973956805, 21.711769035532996],
        [600, 0.519435312355804, 13.930680781758957], [600, 0.5780771438219949, 27.591273615635178], [600, 0.5832500726577585, 32.83358306188925], [600, 0.6271736378897238, 37.97501628664495], [600, 0.6360149795083995, 45.3942996742671],
        [700, 0.5426216540411593, 13.501291571753987], [700, 0.7337019133785521, 23.370683371298405], [700, 0.7729448582372259, 24.18990888382688], [700, 0.8197660974157766, 48.54911161731208], [700, 0.85554445107581, 65.4470159453303]]

for item in data:
    x_data.append(item[0])
    y_data.append(item[1])
    z_data.append(item[2])


parameters, covariance = optimize.curve_fit(function, [x_data, y_data], z_data)
print(parameters)

# create surface function model
# setup data points for calculating surface model
model_x_data = np.linspace(min(x_data), max(x_data), 100)
model_y_data = np.linspace(min(y_data), max(y_data), 100)
# create coordinate arrays for vectorized evaluations
X, Y = np.meshgrid(model_x_data, model_y_data)
# calculate Z coordinate array
Z = function(np.array([X, Y]), *parameters)

z = []
for row in Z:
    row.sort()
    z.append(row)

z = np.array(z)
# print(Z)
# print('####################')
# print(z)

import matplotlib.cm
cm = plt.get_cmap("jet")


# setup figure object
fig = plt.figure()
# setup 3d object
ax = plt.axes(projection='3d')
# plot surface
ax.plot_surface(X, Y, Z, cmap=cm, alpha=0.5)
ax.set_title('TCP Variation in Dissociation Rates')
# plot input data
ax.scatter(x_data, y_data, z_data, color='black')
# set plot descriptions
ax.set_xlabel('Temperature (K)')
ax.invert_xaxis()
ax.set_ylabel('Shear Stress (GPa)')
ax.set_zlabel('Dissociation Rate (per ns)')

plt.show()

# fig = plt.figure()
# ax = plt.axes(projection='3d')
# ax.scatter(xdata, Average_Shear_Stress_List_1GPa, Dissociation_Rates_1GPa)
# ax.scatter(xdata, Average_Shear_Stress_List_2GPa, Dissociation_Rates_2GPa)
# ax.scatter(xdata, Average_Shear_Stress_List_3GPa, Dissociation_Rates_3GPa)
# ax.scatter(xdata, Average_Shear_Stress_List_4GPa, Dissociation_Rates_4GPa)
# ax.scatter(xdata, Average_Shear_Stress_List_5GPa, Dissociation_Rates_5GPa)
# ax.invert_xaxis()
# plt.show()

