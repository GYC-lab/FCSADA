from tkinter import *
from tkinter.filedialog import *
import os
import hashlib
import time
import numpy as np
import math
import matplotlib.pyplot as plt


LOG_LINE_NUM = 0

# 计算零升阻力系数
def get_result_LSZXS(input):
    # 检查输入格式是否符合要求
    output = []
    if len(input) != 3:
        print("Error: input is not right")
    else:
        output1 = input[0] * (input[1] / input[2])
        output.append(output1)
    return output


# 计算翼载荷
def get_result_YZH(input):
    output = []
    if len(input) != 2:
        print("Error: input is not right")
    else:
        output1 = input[0] / input[1]
        output.append(output1)
    return output


# 计算推重比
def get_result_TZB(input):
    output = []
    if len(input) != 2:
        print("Error: input is not right")
    else:
        output1 = input[1] / input[0]
        output.append(output1)
    return output


# 计算功重比
def get_result_GZB(input):
    output = []
    if len(input) != 2:
        print("Error: input is not right")
    else:
        output1 = input[1] / input[0]
        output2 = output1 * 0.73549875
        output.append(output1)
        output.append(output2)
    return output


# 计算推重比与功重比的转化（Gas-Ferrar模型）
def get_result_Gas_Ferrar(input):
    output = []
    if len(input) != 4:
        print("Error: input is not right")
    else:
        output1 = input[0] * input[1] * 0.883 / (input[2] * 745.7 * (input[3] - 0.117))
        output.append(output1)
    return output

# 2.6.2 功重比→推重比
def get_result_GZB_TZB(input):
    output = []
    if len(input) != 4:
        print("Error: input is not right")
    else:
        output1 = input[0] * input[2] * 745.7 * (input[3] - 0.117) / (0.883 * input[1])
        output.append(output1)
    return output

# 3.1.1 机翼展弦比
def get_result_JYZXB(input):
    if len(input) != 2:
        print("Error: input is not right")
    else:
        y = input[0] ** 2 / input[1]
        output = []
        output.append(y)
        return output


# 3.1.2 机翼根梢比
def get_result_JYGSB(input):
    if len(input) != 2:
        print("Error: input is not right")
    else:
        y = input[0] / input[1]
        output = []
        output.append(y)
        return output


# 3.1.3 机翼根弦长和梢弦长
def get_result_GXCHSXC(input):
    if len(input) != 3:
        print("Error: input is not right")
    else:
        y1 = 2 * input[2] / (input[1] + 1) / input[0]
        output = []
        output.append(y1)
        y2 = input[1] * y1
        output.append(y2)
        y3 = 0.5 * input[2] * (y1 + y2)
        output.append(y3)
        return output


# 3.1.4 平均气动弦长
def get_result_PJQDXC(input):
    if len(input) != 2:
        print("Error: input is not right")
    else:
        y = 2 * (input[0] ** 2 + input[0] * input[1] + input[1] ** 2) / 3 / (input[0] + input[1])
        output = []
        output.append(y)
        return output


# 3.1.5 机翼升力线斜率
def get_result_SLXXL(input):
    if len(input) != 4:
        print("Error: input is not right")
    else:
        y = 2 * input[0] * np.pi / (input[0] * (0.5 * (1 / np.cos(input[2])
                                                       + 1 / np.cos(input[3]) + 2 / input[0] / (input[1]) + 1) + 2))
        y2 = 2 * input[0] / (input[0] * (0.5 * (1 / np.cos(input[2])
                                                + 1 / np.cos(input[3]) + 2 / input[0] / (input[1]) + 1) + 2))
        output = []
        output.append(y)
        output.append(y2)
        return output


# 3.2.1 机身长细比
def get_result_JSCXB(input):
    if len(input) != 2:
        print("Error: input is not right")
    else:
        y = input[0] / input[1]
        output = []
        output.append(y)
        return output


# 3.3.1 纵向机身容量参数
def get_result_ZXJSRLCS(input):
    if len(input) != 4:
        print("Error: input is not right")
    else:
        y = input[0] ** 2 * input[1] / input[2] / input[3]
        output = []
        output.append(y)
        return output


# 3.3.2 平尾尾容量
def get_result_PWWRL(input):
    if len(input) != 4:
        print("Error: input is not right")
    else:
        y = input[0] * input[1] / input[2] / input[3]
        output = []
        output.append(y)
        return output


# 3.3.3 平尾展弦比
def get_result_PWZXB(input):
    if len(input) != 2:
        print("Error: input is not right")
    else:
        y = input[0] ** 2 / input[1]
        output = []
        output.append(y)
        return output


# 3.3.4 平尾根梢比
def get_result_PWGSB(input):
    if len(input) != 2:
        print("Error: input is not right")
    else:
        y = input[0] / input[1]
        output = []
        output.append(y)
        return output


# 3.3.5 平尾根弦长和梢弦长
def get_result_PWGXCHSXC(input):
    if len(input) != 3:
        print("Error: input is not right")
    else:
        y1 = 2 * input[2] / (input[1] + 1) / input[0]
        output = []
        output.append(y1)
        y2 = input[1] * y1
        output.append(y2)
        y3 = 0.5 * input[2] * (y1 + y2)
        output.append(y3)
        return output


# 3.4.1 航向机身容量参数
def get_result_HXJSRLCS(input):
    if len(input) != 4:
        print("Error: input is not right")
    else:
        y = input[0] ** 2 * input[1] / input[2] / input[3]
        output = []
        output.append(y)
        return output


# 3.4.2 垂尾尾容量
def get_result_CWWRL(input):
    if len(input) != 4:
        print("Error: input is not right")
    else:
        y = input[0] * input[1] / input[2] / input[3]
        output = []
        output.append(y)
        return output


# 3.4.3 垂尾展弦比
def get_result_CWZXB(input):
    if len(input) != 2:
        print("Error: input is not right")
    else:
        y = input[0] ** 2 / input[1] * 2
        output = []
        output.append(y)
        return output


# 3.4.4 垂尾根梢比
def get_result_CWGSB(input):
    if len(input) != 2:
        print("Error: input is not right")
    else:
        y = input[0] / input[1]
        output = []
        output.append(y)
        return output


# 3.4.5 垂尾根弦长和梢弦长
def get_result_CWGXCHSXC(input):
    if len(input) != 3:
        print("Error: input is not right")
    else:
        y1 = 4 * input[2] / (input[1] + 1) / input[0]
        output = []
        output.append(y1)
        y2 = input[1] * y1
        output.append(y2)
        y3 = 0.5 * input[2] * (y1 + y2)
        output.append(y3)
        return output


# 4.1.1 机翼
def get_result_JYZL(input):
    if len(input) != 6:
        print("Error: input is not right")
    else:
        a = math.cos(input[2] / 180 * np.pi)
        y = 0.00667 * input[0] * (input[1] / a) ** 0.75 \
            * (1 + (1.905 * a / input[1]) ** 0.5) * (input[3] ** 0.55) * (
                        input[1] * input[4] / a / input[0] / input[5]) ** 0.3
        output = []
        output.append(y)
        return output


# 4.1.2 机身
def get_result_JSZL(input):
    if len(input) != 6:
        print("Error: input is not right")
    else:
        y = 1072.6 * ((input[0] * input[1] / 100000) ** 0.286 * (input[2] / 10) ** 0.875 *
                      (input[3] * input[4] / 10) * (input[5] / 100) ** 0.338) ** 1.1
        output = []
        output.append(y)
        return output


# 4.1.3 平尾
def get_result_PW(input):
    if len(input) != 3:
        print("Error: input is not right")
    else:
        y1 = 7.2 * input[0] ** 1.2 * (0.4 + (input[2] + 113) / 935)
        output = []
        output.append(y1)
        y2 = 6.8 * input[1] ** 1.2 * (0.4 + (input[2] + 113) / 1100)
        output.append(y2)
        y3 = y1 + y2
        output.append(y3)
        return output


# 4.1.4 起落架
def get_result_QLJ(input):
    if len(input) != 2:
        print("Error: input is not right")
    else:
        y1 = input[1] * (9.072 + 0.082 * input[0] ** 0.75 + 2.97 * 10 ** (-6) * input[0] ** 1.5)
        output = []
        output.append(y1)
        y2 = input[1] * (18.144 + 0.131 * input[0] ** 0.75 + 0.019 * input[0] + 2.227 * 10 ** (-5) * input[0] ** 1.5)
        output.append(y2)
        y3 = y1 + y2
        output.append(y3)
        return output

#4.1.1 俯冲拉起
def get_result_FCLQ(input):
    if len(input) != 3:
        print("Error: input is not right")
    else:
        y=input[0]**2/input[1]/input[2]+1
        output = []
        output.append(y)
        return output

#4.1.2 等速盘旋
def get_result_DSPX(input):
    if len(input) != 1:
        print("Error: input is not right")
    else:
        y=1/np.cos(input[0]/180*np.pi)
        output = []
        output.append(y)
        return output

#4.1.3 垂直突风
def get_result_CZTF(input):
    if len(input) != 7:
        print("Error: input is not right")
    else:
        mu_g=2*input[0]/input[4]/input[3]/input[2]/input[1]
        K_g=0.88*mu_g/(5.3+mu_g)
        y=1+0.88*K_g*input[2]*input[4]*input[5]*input[6]/2/input[0]
        output = []
        output.append(K_g)
        output.append(y)
        return output

# 翼梁稳定性
def get_result_QD(input):
    if len(input) != 3:
        print("Error: input is not right")
    else:
        sigema=input[0]*input[1]/input[2]
        output = []
        output.append(sigema)
        return output

# 翼梁稳定性
def get_result_WDX(input):
    if len(input) != 6:
        print("Error: input is not right")
    else:
        tau=input[0]/input[1]/input[2]
        tau_cr=(5.35+4/((input[1]/input[3])**2))*(np.pi**2)*input[4]/12/(1-input[5]**2)*((input[2]/input[3])**2)
        output = []
        output.append(tau)
        output.append(tau_cr)
        return output

# 1/n 弦线的后掠角
def get_result_HLJ(input):
    if len(input) != 4:
        print("Error: input is not right")
    else:
        xi = math.atan(math.tan(input[0]/180*np.pi)-4/input[1]/input[2]*(input[3]-1)/(input[3]+1))*180/np.pi
        output = []
        output.append(xi)
        return output

# 机翼升力线斜率
def get_result_SLXXL(input):
    if len(input) != 4:
        print("Error: input is not right")
    else:
        p = 0.5*(1/np.cos(input[1]/180*np.pi)+1/np.cos(input[2]/180*np.pi))+2/input[0]/(input[3]+1)
        CLa=2*np.pi*input[0]/(2+input[0]*p)
        output = []
        output.append(CLa)
        output.append(CLa/np.pi)
        return output

# 奥斯瓦尔德效率因子
def get_result_Oswald(input):
    if len(input) != 1:
        print("Error: input is not right")
    else:
        e = 1.78 *(1-0.045*input[0]**0.68)-0.46
        output = []
        output.append(e)
        return output

# 巡航升力系数
def get_result_PF(input):
    if len(input) != 3:
        print("Error: input is not right")
    else:
        CL=2*input[0]/input[1]/input[2]**2*9.8
        output = []
        output.append(CL)
        return output

# 诱导阻力系数
def get_result_DI(input):
    if len(input) != 3:
        print("Error: input is not right")
    else:
        CD=input[0]**2/np.pi/input[1]/input[2]
        output = []
        output.append(CD)
        return output

# 最大升阻比
def get_result_E_max(input):
    if len(input) != 3:
        print("Error: input is not right")
    else:
        Em=1/2/np.sqrt(input[0]/np.pi/input[1]/input[2])
        output = []
        output.append(Em)
        return output

# 起飞滑跑距离
def get_result_TAKEOFF_distance(input):
    if len(input) != 10:
        print("Error: input is not right")
    else:
        V1=1.2*np.sqrt(2*input[0]/input[2]/input[3]/input[4])
        F0=4*input[5]-input[6]*input[0]
        D=0.5*input[3]*V1**2*input[2]*(input[7]+input[4]**2/np.pi/input[8]/input[9])
        F1=4*input[5]-D
        s=input[0]/2/input[1]*V1**2/(F0-F1)*np.log(F0/F1)
        output = []
        output.append(s)
        return output

# 起飞离地速度
def get_result_TAKEOFF_velocity(input):
    if len(input) != 3:
        print("Error: input is not right")
    else:
        V1=1.2*np.sqrt(2*input[0]*9.8/input[1]/input[2])
        output = []
        output.append(V1)
        return output

# 进场速度
def get_result_LAND_velocity(input):
    if len(input) != 3:
        print("Error: input is not right")
    else:
        V1=1.3*np.sqrt(2*input[0]*9.8/input[1]/input[2])
        output = []
        output.append(V1)
        return output

# 等速巡航航程
def get_result_Cruise_distance(input):
    if len(input) != 5:
        print("Error: input is not right")
    else:
        R=input[0]/input[1]*input[2]*np.log(input[3]/input[4])
        output = []
        output.append(R)
        return output

# 航时
def get_result_Cruise_time(input):
    if len(input) != 4:
        print("Error: input is not right")
    else:
        T=1/input[0]*input[1]*np.log(input[2]/input[3])
        output = []
        output.append(T)
        return output

# 盘旋过载
def get_result_hover(input):
    if len(input) != 3:
        print("Error: input is not right")
    else:
        n=1/np.cos(input[1]/180*np.pi)
        w = input[0]*np.sqrt(n**2-1)/input[2]
        R=input[2]**2/input[0]/np.sqrt(n**2-1)
        t_2pi=2*np.pi/w
        output = []
        output.append(n)
        output.append(w)
        output.append(R)
        output.append(t_2pi)
        return output

# 纵向静稳定裕度
def get_result_H_n(input):
    if len(input) != 3:
        print("Error: input is not right")
    else:
        H_n=(input[0]-input[1])/input[2]
        output = []
        output.append(H_n)
        return output

# 大气参数
def get_result_air(input):
    if len(input) != 3:
        print("Error: input is not right")
    else:
        R=287
        gamma=1.4
        mu_0=1.79e-5
        C=110
        if input[0] < 11000:
            T=288.15-0.0065*input[0]
            p=101325*(1-2.25577*input[0]*1e-5)**5.2588
        else:
            T=216.65
            p=22631*np.exp(-(input[0]-11000)*9.8/287/216.65)
        rho = p/R/T
        q=0.5*rho*input[1]**2
        a=(gamma*R*T)**0.5
        mu=mu_0*(T/288)**1.5*(288+C)/(T+C)
        Ma=input[1]/a
        Re=rho*input[1]*input[2]/mu
        output = []
        output.append(T)
        output.append(T-273.15)
        output.append(rho)
        output.append(p)
        output.append(q)
        output.append(a)
        output.append(mu)
        output.append(Ma)
        output.append(Re)
        output.append(rho/1.225)
        return output

# 翼型升力
def get_result_lift_airfoil(input):
    if len(input) != 3:
        print("Error: input is not right")
    else:
        airfoil_lift=input[0]*input[1]*input[2]
        output = []
        output.append(airfoil_lift)
        return output

# 马赫角
def get_result_mahe_angle(input):
    if len(input) != 1:
        print("Error: input is not right")
    else:
        mu_angle=math.asin(1/(input[0]))*180/np.pi
        output = []
        output.append(mu_angle)
        return output

# 伯努利方程
def get_result_bonuli(input):
    if len(input) != 3:
        print("Error: input is not right")
    else:
        V_inf=(2*(input[0]-input[1])/input[2])**0.5
        output = []
        output.append(V_inf)
        return output

# 线化超音速理论
def get_result_linear(input):
    if len(input) != 2:
        print("Error: input is not right")
    else:
        beta=(input[1]**2-1)**0.5
        c_l=4*input[0]/180*np.pi/beta
        c_d=c_l*input[0]/180*np.pi
        output = []
        output.append(c_l)
        output.append(c_d)
        return output

# 平板层流边界层厚度
def get_result_boundary_layer(input):
    if len(input) != 4:
        print("Error: input is not right")
    else:
        delta=5.48*input[0]/(input[0]*input[1]*input[3]/input[2])**0.5
        output = []
        output.append(delta)
        return output

class MY_GUI():
    def __init__(self, root):
        self.root = root

    def callback_main_qifei(self):
        self.destory_now()
        self.choice = 0
        self.create_now(["空机重量", "乘员重量", "装载重量", "燃油重量"],
                        ["kg", "kg", "kg", "kg"],
                        ["起飞重量"],
                        ["kg"])

    def callback_main_linsheng(self):
        self.destory_now()
        self.choice = 1
        self.create_now(["当量蒙皮摩擦阻力系数", "浸湿面积", "参考面积"],
                        [" ", "m²", "m²"],
                        ["零升阻力系数"],
                        [" "],
                        67, 75)

    def callback_main_yizai(self):
        self.destory_now()
        self.choice = 2
        self.create_now(["飞机重量", "机翼面积"],
                        ["kg", "m²"],
                        ["翼载荷"],
                        ["kg/m²"])

    def callback_main_tuizhong(self):
        self.destory_now()
        self.choice = 3
        self.create_now(["飞机重量", "发动机总推力"],
                        ["kg", "N"],
                        ["推重比"],
                        ["N/kg"])

    def callback_main_gongzhong(self):
        self.destory_now()
        self.choice = 4
        self.create_now(["飞机重量", "发动机总功率"],
                        ["kg", "kW"],
                        ["功重比(kw/kg)","功重比(hp/kg)"],
                        ["kW/kg",'hp/kg'])

    def callback_air_jiyi_zhanxian(self):
        self.destory_now()
        self.choice = 5
        self.create_now(["机翼展长", "机翼面积"],
                        ["m", "m²"],
                        ["机翼展弦比"],
                        [" "])

    def callback_air_jiyi_genshao(self):
        self.destory_now()
        self.choice = 6
        self.create_now(["机翼根弦长", "机翼梢弦长"],
                        ["m", "m"],
                        ["机翼根梢比"],
                        [" "])

    def callback_air_jiyi_xianchang(self):
        self.destory_now()
        self.choice = 7
        self.create_now(["机翼展弦比", "机翼根梢比", "机翼展长"],
                        [" ", " ", "m"],
                        ["机翼梢弦长", "机翼根弦长", "机翼面积"],
                        ["m", "m", "m²"])

    def callback_air_jiyi_qidong(self):
        self.destory_now()
        self.choice = 8
        self.create_now(["机翼根弦长", "机翼梢弦长"],
                        ["m", "m"],
                        ["机翼平均气动弦长"],
                        ["m"])

    def callback_air_jiyi_shengli(self):
        self.destory_now()
        self.choice = 9
        self.create_now(["机翼展弦比", "机翼根梢比", "机翼前缘后掠角", "机翼后缘前掠角"],
                        [" ", " ", "°", "°"],
                        ["机翼升力线斜率", "机翼升力线斜率（*π）"],
                        ["/rad", "π/rad"])

    def callback_air_jishen_changxi(self):
        self.destory_now()
        self.choice = 10
        self.create_now(["机身长度", "机身宽度"],
                        ["m", "m"],
                        ["机身长细比"],
                        [" "])

    def callback_air_pingwei_zongxiang(self):
        self.destory_now()
        self.choice = 11
        self.create_now(["最大机身宽度", "机身长度", "机翼面积", "平均气动弦长"],
                        ["m", "m", "m²", "m"],
                        ["纵向机身容量参数"],
                        [" "])

    def callback_air_pingwei_rongliang(self):
        self.destory_now()
        self.choice = 12
        self.create_now(["平尾尾力臂", "平尾面积", "机翼平均气动弦长", "机翼面积"],
                        ["m", "m²", "m", "m²"],
                        ["平尾尾容量"],
                        [" "])

    def callback_air_pingwei_zhanxian(self):
        self.destory_now()
        self.choice = 13
        self.create_now(["平尾展长", "平尾面积"],
                        ["m", "m²"],
                        ["平尾展弦比"],
                        [" "])

    def callback_air_pingwei_gengshao(self):
        self.destory_now()
        self.choice = 14
        self.create_now(["平尾根弦长", "平尾梢弦长"],
                        ["m", "m"],
                        ["平尾根梢比"],
                        [" "])

    def callback_air_pingwei_xianchang(self):
        self.destory_now()
        self.choice = 15
        self.create_now(["平尾展弦比", "平尾根梢比", "平尾展长"],
                        [" ", " ", "m"],
                        ["平尾梢弦长", "平尾根弦长", "平尾面积"],
                        ["m", "m", "m²"])

    def callback_air_chuiwei_canshu(self):
        self.destory_now()
        self.choice = 16
        self.create_now(["最大机身高度", "机身长度", "机翼面积", "机翼展长"],
                        ["m", "m", "m²", "m"],
                        ["航向机身容量参数"],
                        [" "])

    def callback_air_chuiwei_rongliang(self):
        self.destory_now()
        self.choice = 17
        self.create_now(["垂尾尾力臂", "垂尾面积", "机翼平均气动弦长", "机翼面积"],
                        ["m", "m²", "m", "m²"],
                        ["垂尾尾容量"],
                        [" "])

    def callback_air_chuiwei_zhanxian(self):
        self.destory_now()
        self.choice = 18
        self.create_now(["垂尾(半)展长", "垂尾面积"],
                        ["m", "m²"],
                        ["垂尾展弦比"],
                        [" "])

    def callback_air_chuiwei_gengshao(self):
        self.destory_now()
        self.choice = 19
        self.create_now(["垂尾根弦长", "垂尾梢弦长"],
                        ["m", "m"],
                        ["垂尾根梢比"],
                        [" "])

    def callback_air_chuiwei_xianchang(self):
        self.destory_now()
        self.choice = 20
        self.create_now(["垂尾展弦比", "垂尾根梢比", "垂尾(半)展长"],
                        [" ", " ", "m"],
                        ["垂尾梢弦长", "垂尾根弦长", "垂尾面积"],
                        ["m", "m", "m²"])

    def callback_weight_wing(self):
        self.destory_now()
        self.choice = 21
        self.create_now(["最大零油重量", "翼展", "机翼50%弦线后掠角", "设计过载", "机翼面积", "机翼根部最大厚度"],
                        ["kg", "m", "°", " ", "m²", "m"],
                        ["机翼重量"],
                        ["kg"])

    def callback_weight_fuselage(self):
        self.destory_now()
        self.choice = 22
        self.create_now(["最大起飞重量", "设计过载", "机身长度", "机身最大宽度", "机身最大高度", "设计巡航速度"],
                        ["kg", " ", "m", "m", "m", "km/h"],
                        ["机身重量"],
                        ["kg"])

    def callback_weight_tail(self):
        self.destory_now()
        self.choice = 23
        self.create_now(["平尾面积", "垂尾面积", "设计巡航速度"],
                        ["m²", "m²", 'km/h'],
                        ["平尾重量", '垂尾重量', '尾翼重量'],
                        ["kg", 'kg', 'kg'])

    def callback_weight_landing(self):
        self.destory_now()
        self.choice = 24
        self.create_now(['最大起飞重量', '系数[上单翼-1.0,下单翼-1.08]'],
                        ["kg", " "],
                        ['前起落架重量', '主起落架重量', '起落架重量'],
                        ['kg', 'kg', 'kg'])

    def callback_FCLQ(self):
        self.destory_now()
        self.choice = 25
        self.create_now(['飞行速度', '重力加速度', '俯冲半径'],
                        ["m/s", "m/s²", 'm'],
                        ['过载系数'],
                        [' '])

    def callback_DSPX(self):
        self.destory_now()
        self.choice = 26
        self.create_now(['盘旋倾斜角'],
                        ["°"],
                        ['过载系数'],
                        [' '])

    def callback_CZTF(self):
        self.destory_now()
        self.choice = 27
        self.create_now(['翼载荷', '平均气动弦长', '升力线斜率', '重力加速度', '飞行高度上的空气密度', '飞行速度', '垂直突风速度'],
                        ["N/m²", "m",'/rad','m/s²','kg/m³','m/s','m/s'],
                        ['离散阵风减缓因子','过载系数'],
                        [' ',' '])

    def callback_WDX(self):
        self.destory_now()
        self.choice = 28
        self.create_now(['剪力', '上下缘条形心间的距离', '腹板厚度', '支柱的间距', '弹性模量', '泊松比'],
                        ["N", "m",'m','m²','N/m²',''],
                        ['名义剪应力','临界失稳剪应力'],
                        ['N/m²','N/m²'])

    def callback_QD(self):
        self.destory_now()
        self.choice = 29
        self.create_now(['弯矩', '截面最大高度', '截面惯性矩'],
                        ["N*m", "m",'kg*m²'],
                        ['最大正应力'],
                        ['N/m²'])

    def callback_gongzhong2tuizhong(self):
        self.destory_now()
        self.choice = 30
        self.create_now(['功重比','速度','发动机效率','高度修正因子'],
                        ["hp/kg",'m/s','',''],
                        ['推重比'],
                        ['N/kg'])

    def callback_tuizhong2gongzhong(self):
        self.destory_now()
        self.choice = 31
        self.create_now(['推重比','速度','发动机效率','高度修正因子'],
                        ["N/kg",'m/s','',''],
                        ['功重比'],
                        ['hp/kg'])

    def callback_HLJ(self):
        self.destory_now()
        self.choice = 32
        self.create_now(['前缘后掠角', 'n','展弦比', '根梢比'],
                        ["°", " ",'',''],
                        ['1/n弦线的后掠角'],
                        ['°'])

    def callback_SLXXL(self):
        self.destory_now()
        self.choice = 33
        self.create_now(['机翼展弦比', '机翼前缘后掠角','机翼后缘后掠角', '机翼根梢比'],
                        ["", "°",'°',''],
                        ['机翼升力线斜率','机翼升力线斜率(*pi)'],
                        ['/rad','/rad'])

    def callback_Oswald(self):
        self.destory_now()
        self.choice = 34
        self.create_now(['机翼展弦比'],
                        [""],
                        ['奥斯瓦尔德效率因子'],
                        [''])

    def callback_PF(self):
        self.destory_now()
        self.choice = 35
        self.create_now(['翼载荷','当地空气密度','巡航速度'],
                        ["kg/m²",'kg/m³','m/s'],
                        ['巡航升力系数'],
                        [''])

    def callback_DI(self):
        self.destory_now()
        self.choice = 36
        self.create_now(['升力系数','奥斯瓦尔德效率因子','机翼展弦比'],
                        [" ",'',''],
                        ['诱导阻力系数'],
                        [''])

    def callback_TAKEOFF_distance(self):
        self.destory_now()
        self.choice = 37
        self.create_now(['起飞重量(N)', '重力加速度', '机翼面积','当地空气密度','最大升力系数','发动机推力','地面摩擦阻力系数','零升阻力系数','奥斯瓦尔德效率因子','机翼展弦比'],
                        ["N", 'm/s²', 'm²','kg/m³','','N','','','',''],
                        ['起飞滑跑距离'],
                        ['m'])

    def callback_TAKEOFF_velocity(self):
        self.destory_now()
        self.choice = 38
        self.create_now(['翼载荷','当地空气密度','最大升力系数'],
                        ["kg/m²",'kg/m³',''],
                        ['起飞离地速度'],
                        ['m/s'])

    def callback_LAND_velocity(self):
        self.destory_now()
        self.choice = 39
        self.create_now(['翼载荷','当地空气密度','最大升力系数'],
                        ["kg/m²",'kg/m³',''],
                        ['起飞离地速度'],
                        ['m/s'])

    def callback_Cruise_distance(self):
        self.destory_now()
        self.choice = 40
        self.create_now(['巡航速度','耗油率','升阻比','巡航段初始重量','巡航段结束重量'],
                        ["m/s",'N/(kw·h)','','N','N'],
                        ['等速巡航航程'],
                        ['m'])

    def callback_Cruise_time(self):
        self.destory_now()
        self.choice = 41
        self.create_now(['耗油率','升阻比','巡航段初始重量','巡航段结束重量'],
                        ['N/(kw·h)','','N','N'],
                        ['航时'],
                        ['h'])

    def callback_hover(self):
        self.destory_now()
        self.choice = 42
        self.create_now(['重力加速度','盘旋倾斜角','飞行速度'],
                        ['m/s²','°','m/s'],
                        ['盘旋过载','盘旋角速率','盘旋半径','盘旋一周所需时间'],
                        ['','rad/s','m','s'])

    def callback_E_max(self):
        self.destory_now()
        self.choice = 43
        self.create_now(['零升阻力系数','奥斯瓦尔德效率因子','展弦比'],
                        ['','',''],
                        ['最大升阻比'],
                        [''])

    def callback_H_n(self):
        self.destory_now()
        self.choice = 44
        self.create_now(['全机焦点位置','重心位置','平均气动弦长'],
                        ['m','m','m'],
                        ['纵向静稳定裕度'],
                        [''])

    def callback_air(self):
        self.destory_now()
        self.choice = 45
        self.create_now(['高度', '速度', '特征尺寸'],
                        ['m', 'm/s', 'm'],
                        ['温度(K)','温度(℃)','密度','气压','动压','音速','粘度','马赫数','雷诺数','密度比'],
                        ['T','℃','kg/m³','Pa','N/m²','m/s','N*s/m²','','',''])

    def callback_lift_airfoil(self):
        self.destory_now()
        self.choice = 46
        self.create_now(['空气密度', '来流速度', '速度环量'],
                        ['kg/m³', 'm/s', 'm²/s'],
                        ['翼型升力'],
                        ['N/m'])

    def callback_bonuli(self):
        self.destory_now()
        self.choice = 47
        self.create_now(['总压(皮托管测得压强)', '当前高度大气静压', '当前高度大气密度'],
                        ['Pa', 'Pa', 'kg/m³'],
                        ['飞行速度'],
                        ['m/s'])

    def callback_mahe_angel(self):
        self.destory_now()
        self.choice = 48
        self.create_now(['马赫数'],
                        [''],
                        ['马赫角'],
                        ['°'])

    def callback_linear(self):
        self.destory_now()
        self.choice = 49
        self.create_now(['迎角','马赫数'],
                        ['°',''],
                        ['升力系数','波阻系数'],
                        ['',''])

    def callback_boundary_layer(self):
        self.destory_now()
        self.choice = 50
        self.create_now(['距离平板前缘距离','气体密度','气体粘度','来流速度'],
                        ['m','kg/m³','N*s/m²','m/s'],
                        ['边界层厚度'],
                        ['m'])

    def callback(self):
        print('switch')

    # 将当前输入输出的控件清除
    def destory_now(self):
        for index in range(0, len(self.inputLabel)):
            if self.inputLabel[index] is None:
                break
            self.inputLabel[index].destroy()
            self.inputLabel[index] = None
            self.inputText[index].destroy()
            self.inputText[index] = None
            self.inputUnit[index].destroy()
            self.inputUnit[index] = None

        for index in range(0, len(self.outputLabel)):
            if self.outputLabel[index] is None:
                break
            self.outputLabel[index].destroy()
            self.outputLabel[index] = None
            self.outputText[index].destroy()
            self.outputText[index] = None
            self.outputUnit[index].destroy()
            self.outputUnit[index] = None

    # 根据输入的四个数组以及两个长度，创建对应模块需求的控件
    # inLabel--输入区域的标签 inUnit--输入区域的单位 outLabel--输出区域的标签 outUnit--输出区域的单位
    # inLength--输入区域的长度 outLength--输出区域的长度
    def create_now(self, inLabel, inUnit, outLabel, outUnit, inLength=75, outLength=75):
        if len(inLabel) != len(inUnit) or len(outLabel) != len(outUnit):
            print("Error: the length of label and unit of input or output not equal")
        for index in range(0, len(inLabel)):
            self.inputLabel[index] = Label(self.inputFrame, text=inLabel[index], font='50px')
            self.inputText[index] = Entry(self.inputFrame, width=inLength, font='50px')
            self.inputUnit[index] = Label(self.inputFrame, text=inUnit[index], font='50px')
        for index in range(0, len(outLabel)):
            self.outputLabel[index] = Label(self.outputFrame, text=outLabel[index], font='50px')
            self.outputText[index] = Text(self.outputFrame, width=outLength, height=1, font='50px', bg='#FFFFFF',
                                          highlightthickness=1,
                                          relief='flat')
            self.outputUnit[index] = Label(self.outputFrame, text=outUnit[index], font='50px')
        self.show_input_lable()
        self.show_input_text()
        self.show_input_unit()
        self.show_output_lable()
        self.show_output_text()
        self.show_output_unit()

    # 绘制输入框标签
    def show_input_lable(self):
        showRow = 1
        for input_label in self.inputLabel:
            if input_label is None:
                break
            input_label.grid(row=showRow, column=1, pady=5)
            showRow += 1

    # 绘制输入框
    def show_input_text(self):
        showRow = 1
        for input_text in self.inputText:
            if input_text is None:
                break
            input_text.grid(row=showRow, column=3, rowspan=1, columnspan=1, padx=10)
            showRow += 1

    # 绘制输入单位
    def show_input_unit(self):
        showRow = 1
        for input_unit in self.inputUnit:
            if input_unit is None:
                break
            input_unit.grid(row=showRow, column=4, rowspan=1, columnspan=1, padx=10)
            showRow += 1

    # 绘制输出框标签
    def show_output_lable(self):
        showRow = 1
        for output_label in self.outputLabel:
            if output_label is None:
                break
            output_label.grid(row=showRow, column=1, pady=5)
            showRow += 1

    # 绘制输入框
    def show_output_text(self):
        showRow = 1
        for output_text in self.outputText:
            if output_text is None:
                break
            output_text.grid(row=showRow, column=3, rowspan=1, columnspan=1, padx=10)
            showRow += 1

    # 绘制输入单位
    def show_output_unit(self):
        showRow = 1
        for output_unit in self.outputUnit:
            if output_unit is None:
                break
            output_unit.grid(row=showRow, column=4, rowspan=1, columnspan=1, padx=10)
            showRow += 1

    # 获取计算结果并将其记录在历史中
    def get_result(self):
        if self.choice == 0:
            try:
                input1 = float(self.inputText[0].get())
            except:
                self.logText.insert(1.0, "错误：第一个参数存在格式问题，请检查\n")
                return
            try:
                input2 = float(self.inputText[1].get())
            except:
                self.logText.insert(1.0, "错误：第二个参数存在格式问题，请检查\n")
                return
            try:
                input3 = float(self.inputText[2].get())
            except:
                self.logText.insert(1.0, "错误：第三个参数存在格式问题，请检查\n")
                return
            try:
                input4 = float(self.inputText[3].get())
            except:
                self.logText.insert(1.0, "错误：第四个参数存在格式问题，请检查\n")
                return

            result1 = input1 + input2 + input3 + input4
            self.outputText[0].delete(1.0, END)
            self.outputText[0].insert(1.0, result1)
            self.logText.insert(1.0, "----------成功----------\n项目：总体设计-起飞重量\n参数为" +
                                "\n\t空机重量: " + str(input1) + "kg" +
                                "\n\t乘员重量: " + str(input2) + "kg" +
                                "\n\t装载重量: " + str(input3) + "kg" +
                                "\n\t燃油重量: " + str(input4) + "kg" +
                                "\n结果为" +
                                "\n\t起飞重量：" + str(result1) + "kg" +
                                "\n")
            self.history.append([0, input1, input2, input3, input4, result1])
            print(len(self.history))
            self.nowHistory = len(self.history) - 1

        elif self.choice == 1:
            try:
                input1 = float(self.inputText[0].get())
            except:
                self.logText.insert(1.0, "错误：第一个参数存在格式问题，请检查\n")
                return
            try:
                input2 = float(self.inputText[1].get())
            except:
                self.logText.insert(1.0, "错误：第二个参数存在格式问题，请检查\n")
                return
            try:
                input3 = float(self.inputText[2].get())
            except:
                self.logText.insert(1.0, "错误：第三个参数存在格式问题，请检查\n")
                return

            result = get_result_LSZXS([input1, input2, input3])

            self.outputText[0].delete(1.0, END)
            self.outputText[0].insert(1.0, result[0])
            self.logText.insert(1.0, "----------成功----------\n项目：总体设计-零升阻力系数\n参数为" +
                                "\n\t当量蒙皮摩擦阻力系数: " + str(input1) + " " +
                                "\n\t浸湿面积: " + str(input2) + "m²" +
                                "\n\t参考面积: " + str(input3) + "m²" +
                                "\n结果为" +
                                "\n\t零升阻力系数：" + str(result[0]) + " " +
                                "\n")
            self.history.append([0, input1, input2, input3, result[0]])
            print(len(self.history))
            self.nowHistory = len(self.history) - 1

        elif self.choice == 2:
            try:
                input1 = float(self.inputText[0].get())
            except:
                self.logText.insert(1.0, "错误：第一个参数存在格式问题，请检查\n")
                return
            try:
                input2 = float(self.inputText[1].get())
            except:
                self.logText.insert(1.0, "错误：第二个参数存在格式问题，请检查\n")
                return

            result = get_result_YZH([input1, input2])

            self.outputText[0].delete(1.0, END)
            self.outputText[0].insert(1.0, result[0])
            self.logText.insert(1.0, "----------成功----------\n项目：总体设计-翼载荷\n参数为" +
                                "\n\t飞机重量: " + str(input1) + "kg" +
                                "\n\t机翼面积: " + str(input2) + "m²" +
                                "\n结果为" +
                                "\n\t翼载荷：" + str(result[0]) + "kg/m²" +
                                "\n")
            self.history.append([0, input1, input2, result[0]])
            print(len(self.history))
            self.nowHistory = len(self.history) - 1

        elif self.choice == 3:
            try:
                input1 = float(self.inputText[0].get())
            except:
                self.logText.insert(1.0, "错误：第一个参数存在格式问题，请检查\n")
                return
            try:
                input2 = float(self.inputText[1].get())
            except:
                self.logText.insert(1.0, "错误：第二个参数存在格式问题，请检查\n")
                return

            result = get_result_TZB([input1, input2])

            self.outputText[0].delete(1.0, END)
            self.outputText[0].insert(1.0, result[0])
            self.logText.insert(1.0, "----------成功----------\n项目：总体设计-推重比\n参数为" +
                                "\n\t飞机重量: " + str(input1) + "kg" +
                                "\n\t发动机总推力: " + str(input2) + "N" +
                                "\n结果为" +
                                "\n\t推重比：" + str(result[0]) + "N/kg" +
                                "\n")
            self.history.append([0, input1, input2, result[0]])
            print(len(self.history))
            self.nowHistory = len(self.history) - 1

        elif self.choice == 4:
            try:
                input1 = float(self.inputText[0].get())
            except:
                self.logText.insert(1.0, "错误：第一个参数存在格式问题，请检查\n")
                return
            try:
                input2 = float(self.inputText[1].get())
            except:
                self.logText.insert(1.0, "错误：第二个参数存在格式问题，请检查\n")
                return

            result = get_result_GZB([input1, input2])

            self.outputText[0].delete(1.0, END)
            self.outputText[0].insert(1.0, result[0])
            self.outputText[1].delete(1.0, END)
            self.outputText[1].insert(1.0, result[1])
            self.logText.insert(1.0, "----------成功----------\n项目：总体设计-功重比\n参数为" +
                                "\n\t飞机重量: " + str(input1) + "kg" +
                                "\n\t发动机总功率: " + str(input2) + "kW" +
                                "\n结果为" +
                                "\n\t功重比(kw/kg)：" + str(result[0]) + "kW/kg" +
                                "\n\t功重比(hp/kg)：" + str(result[1]) + "hp/kg" +
                                "\n")
            self.history.append([0, input1, input2, result[0],result[1]])
            print(len(self.history))
            self.nowHistory = len(self.history) - 1

        elif self.choice == 5:
            try:
                input1 = float(self.inputText[0].get())
            except:
                self.logText.insert(1.0, "错误：第一个参数存在格式问题，请检查\n")
                return
            try:
                input2 = float(self.inputText[1].get())
            except:
                self.logText.insert(1.0, "错误：第二个参数存在格式问题，请检查\n")
                return

            result = get_result_JYZXB([input1, input2])

            self.outputText[0].delete(1.0, END)
            self.outputText[0].insert(1.0, result[0])
            self.logText.insert(1.0, "----------成功----------\n项目：气动设计-展弦比\n参数为" +
                                "\n\t机翼展长: " + str(input1) + "m" +
                                "\n\t机翼面积: " + str(input2) + "m²" +
                                "\n结果为" +
                                "\n\t展弦比：" + str(result[0]) + " " +
                                "\n")
            self.history.append([0, input1, input2, result[0]])
            print(len(self.history))
            self.nowHistory = len(self.history) - 1

        elif self.choice == 6:
            try:
                input1 = float(self.inputText[0].get())
            except:
                self.logText.insert(1.0, "错误：第一个参数存在格式问题，请检查\n")
                return
            try:
                input2 = float(self.inputText[1].get())
            except:
                self.logText.insert(1.0, "错误：第二个参数存在格式问题，请检查\n")
                return

            result = get_result_JYGSB([input1, input2])

            self.outputText[0].delete(1.0, END)
            self.outputText[0].insert(1.0, result[0])
            self.logText.insert(1.0, "----------成功----------\n项目：气动设计-机翼根梢比\n参数为" +
                                "\n\t机翼根弦长: " + str(input1) + "m" +
                                "\n\t机翼梢弦长: " + str(input2) + "m" +
                                "\n结果为" +
                                "\n\t机翼根梢比：" + str(result[0]) + " " +
                                "\n")
            self.history.append([0, input1, input2, result[0]])
            print(len(self.history))
            self.nowHistory = len(self.history) - 1

        elif self.choice == 7:
            try:
                input1 = float(self.inputText[0].get())
            except:
                self.logText.insert(1.0, "错误：第一个参数存在格式问题，请检查\n")
                return
            try:
                input2 = float(self.inputText[1].get())
            except:
                self.logText.insert(1.0, "错误：第二个参数存在格式问题，请检查\n")
                return
            try:
                input3 = float(self.inputText[2].get())
            except:
                self.logText.insert(1.0, "错误：第三个参数存在格式问题，请检查\n")
                return

            result = get_result_GXCHSXC([input1, input2, input3])

            self.outputText[0].delete(1.0, END)
            self.outputText[0].insert(1.0, result[0])
            self.outputText[1].delete(1.0, END)
            self.outputText[1].insert(1.0, result[1])
            self.outputText[2].delete(1.0, END)
            self.outputText[2].insert(1.0, result[2])
            self.logText.insert(1.0, "----------成功----------\n项目：气动设计-机翼根弦长和梢弦长\n参数为" +
                                "\n\t机翼展弦比: " + str(input1) + " " +
                                "\n\t机翼根梢比: " + str(input2) + " " +
                                "\n\t机翼展长: " + str(input3) + "m" +
                                "\n结果为" +
                                "\n\t机翼梢弦长：" + str(result[0]) + "m" +
                                "\n\t机翼根弦长：" + str(result[1]) + "m" +
                                "\n\t机翼面积：" + str(result[2]) + "m²" +
                                "\n")
            self.history.append([0, input1, input2, input3, result[0], result[1], result[2]])
            print(len(self.history))
            self.nowHistory = len(self.history) - 1

        elif self.choice == 8:
            try:
                input1 = float(self.inputText[0].get())
            except:
                self.logText.insert(1.0, "错误：第一个参数存在格式问题，请检查\n")
                return
            try:
                input2 = float(self.inputText[1].get())
            except:
                self.logText.insert(1.0, "错误：第二个参数存在格式问题，请检查\n")
                return

            result = get_result_PJQDXC([input1, input2])

            self.outputText[0].delete(1.0, END)
            self.outputText[0].insert(1.0, result[0])
            self.logText.insert(1.0, "----------成功----------\n项目：气动设计-机翼平均气动弦长\n参数为" +
                                "\n\t机翼根弦长: " + str(input1) + "m" +
                                "\n\t机翼梢弦长: " + str(input2) + "m" +
                                "\n结果为" +
                                "\n\t机翼平均气动弦长：" + str(result[0]) + "m" +
                                "\n")
            self.history.append([0, input1, input2, result[0]])
            print(len(self.history))
            self.nowHistory = len(self.history) - 1

        elif self.choice == 9:
            try:
                input1 = float(self.inputText[0].get())
            except:
                self.logText.insert(1.0, "错误：第一个参数存在格式问题，请检查\n")
                return
            try:
                input2 = float(self.inputText[1].get())
            except:
                self.logText.insert(1.0, "错误：第二个参数存在格式问题，请检查\n")
                return
            try:
                input3 = float(self.inputText[2].get())
            except:
                self.logText.insert(1.0, "错误：第三个参数存在格式问题，请检查\n")
                return
            try:
                input4 = float(self.inputText[3].get())
            except:
                self.logText.insert(1.0, "错误：第四个参数存在格式问题，请检查\n")
                return

            result = get_result_SLXXL([input1, input2, input3, input4])

            self.outputText[0].delete(1.0, END)
            self.outputText[0].insert(1.0, result[0])
            self.outputText[1].delete(1.0, END)
            self.outputText[1].insert(1.0, result[1])
            self.logText.insert(1.0, "----------成功----------\n项目：气动设计-机翼升力线斜率\n参数为" +
                                "\n\t机翼展弦比: " + str(input1) + "" +
                                "\n\t机翼根梢比: " + str(input2) + "" +
                                "\n\t机翼前缘后掠角: " + str(input3) + "" +
                                "\n\t机翼后缘前掠角: " + str(input4) + "" +
                                "\n结果为" +
                                "\n\t机翼升力线斜率：" + str(result[0]) + "/rad" +
                                "\n\t机翼升力线斜率(*π)：" + str(result[1]) + "π/rad" +
                                "\n")
            self.history.append([0, input1, input2, input3, input4, result[0], result[1]])
            print(len(self.history))
            self.nowHistory = len(self.history) - 1

        elif self.choice == 10:
            try:
                input1 = float(self.inputText[0].get())
            except:
                self.logText.insert(1.0, "错误：第一个参数存在格式问题，请检查\n")
                return
            try:
                input2 = float(self.inputText[1].get())
            except:
                self.logText.insert(1.0, "错误：第二个参数存在格式问题，请检查\n")
                return

            result = get_result_JSCXB([input1, input2])

            self.outputText[0].delete(1.0, END)
            self.outputText[0].insert(1.0, result[0])
            self.logText.insert(1.0, "----------成功----------\n项目：气动设计-机身长细比\n参数为" +
                                "\n\t机身长度: " + str(input1) + "m" +
                                "\n\t机身高度: " + str(input2) + "m" +
                                "\n结果为" +
                                "\n\t机身长细比：" + str(result[0]) + "" +
                                "\n")
            self.history.append([0, input1, input2, result[0]])
            print(len(self.history))
            self.nowHistory = len(self.history) - 1

        elif self.choice == 11:
            try:
                input1 = float(self.inputText[0].get())
            except:
                self.logText.insert(1.0, "错误：第一个参数存在格式问题，请检查\n")
                return
            try:
                input2 = float(self.inputText[1].get())
            except:
                self.logText.insert(1.0, "错误：第二个参数存在格式问题，请检查\n")
                return
            try:
                input3 = float(self.inputText[2].get())
            except:
                self.logText.insert(1.0, "错误：第三个参数存在格式问题，请检查\n")
                return
            try:
                input4 = float(self.inputText[3].get())
            except:
                self.logText.insert(1.0, "错误：第四个参数存在格式问题，请检查\n")
                return

            result = get_result_ZXJSRLCS([input1, input2, input3, input4])

            self.outputText[0].delete(1.0, END)
            self.outputText[0].insert(1.0, result[0])
            self.logText.insert(1.0, "----------成功----------\n项目：气动设计-纵向机身容量参数\n参数为" +
                                "\n\t最大机身宽度: " + str(input1) + "m" +
                                "\n\t机身长度: " + str(input2) + "m" +
                                "\n\t机翼面积: " + str(input3) + "m²" +
                                "\n\t平均气动弦长: " + str(input4) + "m" +
                                "\n结果为" +
                                "\n\t纵向机身容量参数：" + str(result[0]) + "" +
                                "\n")
            self.history.append([0, input1, input2, input3, input4, result[0]])
            print(len(self.history))
            self.nowHistory = len(self.history) - 1

        elif self.choice == 12:
            try:
                input1 = float(self.inputText[0].get())
            except:
                self.logText.insert(1.0, "错误：第一个参数存在格式问题，请检查\n")
                return
            try:
                input2 = float(self.inputText[1].get())
            except:
                self.logText.insert(1.0, "错误：第二个参数存在格式问题，请检查\n")
                return
            try:
                input3 = float(self.inputText[2].get())
            except:
                self.logText.insert(1.0, "错误：第三个参数存在格式问题，请检查\n")
                return
            try:
                input4 = float(self.inputText[3].get())
            except:
                self.logText.insert(1.0, "错误：第四个参数存在格式问题，请检查\n")
                return

            result = get_result_PWWRL([input1, input2, input3, input4])

            self.outputText[0].delete(1.0, END)
            self.outputText[0].insert(1.0, result[0])
            self.logText.insert(1.0, "----------成功----------\n项目：气动设计-平尾尾容量\n参数为" +
                                "\n\t平尾尾力臂: " + str(input1) + "m" +
                                "\n\t平尾面积: " + str(input2) + "m²" +
                                "\n\t机翼平均气动弦长: " + str(input3) + "m" +
                                "\n\t机翼面积: " + str(input4) + "m²" +
                                "\n结果为" +
                                "\n\t平尾尾容量：" + str(result[0]) + "" +
                                "\n")
            self.history.append([0, input1, input2, input3, input4, result[0]])
            print(len(self.history))
            self.nowHistory = len(self.history) - 1

        elif self.choice == 13:
            try:
                input1 = float(self.inputText[0].get())
            except:
                self.logText.insert(1.0, "错误：第一个参数存在格式问题，请检查\n")
                return
            try:
                input2 = float(self.inputText[1].get())
            except:
                self.logText.insert(1.0, "错误：第二个参数存在格式问题，请检查\n")
                return

            result = get_result_PWZXB([input1, input2])

            self.outputText[0].delete(1.0, END)
            self.outputText[0].insert(1.0, result[0])
            self.logText.insert(1.0, "----------成功----------\n项目：气动设计-平尾展弦比\n参数为" +
                                "\n\t平尾展长: " + str(input1) + "m" +
                                "\n\t平尾面积: " + str(input2) + "m²" +
                                "\n结果为" +
                                "\n\t平尾展弦比：" + str(result[0]) + "" +
                                "\n")
            self.history.append([0, input1, input2, result[0]])
            print(len(self.history))
            self.nowHistory = len(self.history) - 1

        elif self.choice == 14:
            try:
                input1 = float(self.inputText[0].get())
            except:
                self.logText.insert(1.0, "错误：第一个参数存在格式问题，请检查\n")
                return
            try:
                input2 = float(self.inputText[1].get())
            except:
                self.logText.insert(1.0, "错误：第二个参数存在格式问题，请检查\n")
                return

            result = get_result_PWGSB([input1, input2])

            self.outputText[0].delete(1.0, END)
            self.outputText[0].insert(1.0, result[0])
            self.logText.insert(1.0, "----------成功----------\n项目：气动设计-平尾根梢比\n参数为" +
                                "\n\t平尾根弦长: " + str(input1) + "m" +
                                "\n\t平尾梢弦长: " + str(input2) + "m" +
                                "\n结果为" +
                                "\n\t平尾根梢比：" + str(result[0]) + "" +
                                "\n")
            self.history.append([0, input1, input2, result[0]])
            print(len(self.history))
            self.nowHistory = len(self.history) - 1

        elif self.choice == 15:
            try:
                input1 = float(self.inputText[0].get())
            except:
                self.logText.insert(1.0, "错误：第一个参数存在格式问题，请检查\n")
                return
            try:
                input2 = float(self.inputText[1].get())
            except:
                self.logText.insert(1.0, "错误：第二个参数存在格式问题，请检查\n")
                return
            try:
                input3 = float(self.inputText[2].get())
            except:
                self.logText.insert(1.0, "错误：第三个参数存在格式问题，请检查\n")
                return

            result = get_result_PWGXCHSXC([input1, input2, input3])

            self.outputText[0].delete(1.0, END)
            self.outputText[0].insert(1.0, result[0])
            self.outputText[1].delete(1.0, END)
            self.outputText[1].insert(1.0, result[1])
            self.outputText[2].delete(1.0, END)
            self.outputText[2].insert(1.0, result[2])
            self.logText.insert(1.0, "----------成功----------\n项目：气动设计-平尾根弦长和梢弦长\n参数为" +
                                "\n\t平尾展弦比: " + str(input1) + " " +
                                "\n\t平尾根梢比: " + str(input2) + " " +
                                "\n\t平尾展长: " + str(input3) + "m" +
                                "\n结果为" +
                                "\n\t平尾梢弦长：" + str(result[0]) + "m" +
                                "\n\t平尾根弦长：" + str(result[1]) + "m" +
                                "\n\t平尾面积：" + str(result[2]) + "m²" +
                                "\n")
            self.history.append([0, input1, input2, input3, result[0], result[1], result[2]])
            print(len(self.history))
            self.nowHistory = len(self.history) - 1

        elif self.choice == 16:
            try:
                input1 = float(self.inputText[0].get())
            except:
                self.logText.insert(1.0, "错误：第一个参数存在格式问题，请检查\n")
                return
            try:
                input2 = float(self.inputText[1].get())
            except:
                self.logText.insert(1.0, "错误：第二个参数存在格式问题，请检查\n")
                return
            try:
                input3 = float(self.inputText[2].get())
            except:
                self.logText.insert(1.0, "错误：第三个参数存在格式问题，请检查\n")
                return
            try:
                input4 = float(self.inputText[3].get())
            except:
                self.logText.insert(1.0, "错误：第四个参数存在格式问题，请检查\n")
                return

            result = get_result_HXJSRLCS([input1, input2, input3, input4])

            self.outputText[0].delete(1.0, END)
            self.outputText[0].insert(1.0, result[0])
            self.logText.insert(1.0, "----------成功----------\n项目：气动设计-航向机身容量参数\n参数为" +
                                "\n\t最大机身高度: " + str(input1) + "m" +
                                "\n\t机身长度: " + str(input2) + "m" +
                                "\n\t机翼面积: " + str(input3) + "m²" +
                                "\n\t机翼展长: " + str(input4) + "m" +
                                "\n结果为" +
                                "\n\t航向机身容量参数：" + str(result[0]) + "" +
                                "\n")
            self.history.append([0, input1, input2, input3, input4, result[0]])
            print(len(self.history))
            self.nowHistory = len(self.history) - 1

        elif self.choice == 17:
            try:
                input1 = float(self.inputText[0].get())
            except:
                self.logText.insert(1.0, "错误：第一个参数存在格式问题，请检查\n")
                return
            try:
                input2 = float(self.inputText[1].get())
            except:
                self.logText.insert(1.0, "错误：第二个参数存在格式问题，请检查\n")
                return
            try:
                input3 = float(self.inputText[2].get())
            except:
                self.logText.insert(1.0, "错误：第三个参数存在格式问题，请检查\n")
                return
            try:
                input4 = float(self.inputText[3].get())
            except:
                self.logText.insert(1.0, "错误：第四个参数存在格式问题，请检查\n")
                return

            result = get_result_CWWRL([input1, input2, input3, input4])

            self.outputText[0].delete(1.0, END)
            self.outputText[0].insert(1.0, result[0])
            self.logText.insert(1.0, "----------成功----------\n项目：气动设计-垂尾尾容量\n参数为" +
                                "\n\t垂尾尾力臂: " + str(input1) + "m" +
                                "\n\t垂尾面积: " + str(input2) + "m²" +
                                "\n\t机翼平均气动弦长: " + str(input3) + "m" +
                                "\n\t机翼面积: " + str(input4) + "m²" +
                                "\n结果为" +
                                "\n\t垂尾尾容量：" + str(result[0]) + "" +
                                "\n")
            self.history.append([0, input1, input2, input3, input4, result[0]])
            print(len(self.history))
            self.nowHistory = len(self.history) - 1

        elif self.choice == 18:
            try:
                input1 = float(self.inputText[0].get())
            except:
                self.logText.insert(1.0, "错误：第一个参数存在格式问题，请检查\n")
                return
            try:
                input2 = float(self.inputText[1].get())
            except:
                self.logText.insert(1.0, "错误：第二个参数存在格式问题，请检查\n")
                return

            result = get_result_CWZXB([input1, input2])

            self.outputText[0].delete(1.0, END)
            self.outputText[0].insert(1.0, result[0])
            self.logText.insert(1.0, "----------成功----------\n项目：气动设计-垂尾展弦比\n参数为" +
                                "\n\t垂尾(半)展长: " + str(input1) + "m" +
                                "\n\t垂尾面积: " + str(input2) + "m²" +
                                "\n结果为" +
                                "\n\t垂尾展弦比：" + str(result[0]) + "" +
                                "\n")
            self.history.append([0, input1, input2, result[0]])
            print(len(self.history))
            self.nowHistory = len(self.history) - 1

        elif self.choice == 19:
            try:
                input1 = float(self.inputText[0].get())
            except:
                self.logText.insert(1.0, "错误：第一个参数存在格式问题，请检查\n")
                return
            try:
                input2 = float(self.inputText[1].get())
            except:
                self.logText.insert(1.0, "错误：第二个参数存在格式问题，请检查\n")
                return

            result = get_result_CWGSB([input1, input2])

            self.outputText[0].delete(1.0, END)
            self.outputText[0].insert(1.0, result[0])
            self.logText.insert(1.0, "----------成功----------\n项目：气动设计-垂尾根梢比\n参数为" +
                                "\n\t垂尾根弦长: " + str(input1) + "m" +
                                "\n\t垂尾梢弦长: " + str(input2) + "m" +
                                "\n结果为" +
                                "\n\t垂尾根梢比：" + str(result[0]) + "" +
                                "\n")
            self.history.append([0, input1, input2, result[0]])
            print(len(self.history))
            self.nowHistory = len(self.history) - 1

        elif self.choice == 20:
            try:
                input1 = float(self.inputText[0].get())
            except:
                self.logText.insert(1.0, "错误：第一个参数存在格式问题，请检查\n")
                return
            try:
                input2 = float(self.inputText[1].get())
            except:
                self.logText.insert(1.0, "错误：第二个参数存在格式问题，请检查\n")
                return
            try:
                input3 = float(self.inputText[2].get())
            except:
                self.logText.insert(1.0, "错误：第三个参数存在格式问题，请检查\n")
                return

            result = get_result_CWGXCHSXC([input1, input2, input3])

            self.outputText[0].delete(1.0, END)
            self.outputText[0].insert(1.0, result[0])
            self.outputText[1].delete(1.0, END)
            self.outputText[1].insert(1.0, result[1])
            self.outputText[2].delete(1.0, END)
            self.outputText[2].insert(1.0, result[2])
            self.logText.insert(1.0, "----------成功----------\n项目：气动设计-垂尾根弦长和梢弦长\n参数为" +
                                "\n\t垂尾展弦比: " + str(input1) + " " +
                                "\n\t垂尾根梢比: " + str(input2) + " " +
                                "\n\t垂尾(半)展长: " + str(input3) + "m" +
                                "\n结果为" +
                                "\n\t垂尾梢弦长：" + str(result[0]) + "m" +
                                "\n\t垂尾根弦长：" + str(result[1]) + "m" +
                                "\n\t垂尾面积：" + str(result[2]) + "m²" +
                                "\n")
            self.history.append([0, input1, input2, input3, result[0], result[1], result[2]])
            print(len(self.history))
            self.nowHistory = len(self.history) - 1

        elif self.choice == 21:
            try:
                input1 = float(self.inputText[0].get())
            except:
                self.logText.insert(1.0, "错误：第一个参数存在格式问题，请检查\n")
                return
            try:
                input2 = float(self.inputText[1].get())
            except:
                self.logText.insert(1.0, "错误：第二个参数存在格式问题，请检查\n")
                return
            try:
                input3 = float(self.inputText[2].get())
            except:
                self.logText.insert(1.0, "错误：第三个参数存在格式问题，请检查\n")
                return
            try:
                input4 = float(self.inputText[3].get())
            except:
                self.logText.insert(1.0, "错误：第四个参数存在格式问题，请检查\n")
                return
            try:
                input5 = float(self.inputText[4].get())
            except:
                self.logText.insert(1.0, "错误：第五个参数存在格式问题，请检查\n")
                return
            try:
                input6 = float(self.inputText[5].get())
            except:
                self.logText.insert(1.0, "错误：第六个参数存在格式问题，请检查\n")
                return

            result = get_result_JYZL([input1, input2, input3, input4, input5, input6])

            self.outputText[0].delete(1.0, END)
            self.outputText[0].insert(1.0, result[0])
            self.logText.insert(1.0, "----------成功----------\n项目：重量分析-机翼重量估算\n参数为" +
                                "\n\t最大零油重量: " + str(input1) + "kg" +
                                "\n\t翼展: " + str(input2) + "m" +
                                "\n\t机翼50%弦线后掠角: " + str(input3) + "°" +
                                "\n\t设计过载: " + str(input4) + " " +
                                "\n\t机翼面积: " + str(input5) + "m²" +
                                "\n\t机翼根部最大厚度: " + str(input6) + "m" +
                                "\n结果为" +
                                "\n\t机翼重量：" + str(result[0]) + "kg" +
                                "\n")
            self.history.append([0, input1, input2, input3, input4, input5, input6, result[0]])
            print(len(self.history))
            self.nowHistory = len(self.history) - 1

        elif self.choice == 22:
            try:
                input1 = float(self.inputText[0].get())
            except:
                self.logText.insert(1.0, "错误：第一个参数存在格式问题，请检查\n")
                return
            try:
                input2 = float(self.inputText[1].get())
            except:
                self.logText.insert(1.0, "错误：第二个参数存在格式问题，请检查\n")
                return
            try:
                input3 = float(self.inputText[2].get())
            except:
                self.logText.insert(1.0, "错误：第三个参数存在格式问题，请检查\n")
                return
            try:
                input4 = float(self.inputText[3].get())
            except:
                self.logText.insert(1.0, "错误：第四个参数存在格式问题，请检查\n")
                return
            try:
                input5 = float(self.inputText[4].get())
            except:
                self.logText.insert(1.0, "错误：第五个参数存在格式问题，请检查\n")
                return
            try:
                input6 = float(self.inputText[5].get())
            except:
                self.logText.insert(1.0, "错误：第六个参数存在格式问题，请检查\n")
                return

            result = get_result_JSZL([input1, input2, input3, input4, input5, input6])

            self.outputText[0].delete(1.0, END)
            self.outputText[0].insert(1.0, result[0])
            self.logText.insert(1.0, "----------成功----------\n项目：重量分析-机翼重量估算\n参数为" +
                                "\n\t最大机身重量: " + str(input1) + "kg" +
                                "\n\t设计过载: " + str(input2) + "" +
                                "\n\t机身长度: " + str(input3) + "m" +
                                "\n\t机身最大宽度: " + str(input4) + "m" +
                                "\n\t机身最大高度: " + str(input5) + "m" +
                                "\n\t设计巡航速度: " + str(input6) + "m/s" +
                                "\n结果为" +
                                "\n\t机身重量：" + str(result[0]) + "kg" +
                                "\n")
            self.history.append([0, input1, input2, input3, input4, input5, input6, result[0]])
            print(len(self.history))
            self.nowHistory = len(self.history) - 1

        elif self.choice == 24:
            try:
                input1 = float(self.inputText[0].get())
            except:
                self.logText.insert(1.0, "错误：第一个参数存在格式问题，请检查\n")
                return
            try:
                input2 = float(self.inputText[1].get())
            except:
                self.logText.insert(1.0, "错误：第二个参数存在格式问题，请检查\n")
                return

            result = get_result_QLJ([input1, input2])

            self.outputText[0].delete(1.0, END)
            self.outputText[0].insert(1.0, result[0])
            self.outputText[1].delete(1.0, END)
            self.outputText[1].insert(1.0, result[1])
            self.outputText[2].delete(1.0, END)
            self.outputText[2].insert(1.0, result[2])
            self.logText.insert(1.0, "----------成功----------\n项目：重量分析-起落架重量估算\n参数为" +
                                "\n\t最大起飞重量: " + str(input1) + "kg" +
                                "\n\t比例系数: " + str(input2) + "" +
                                "\n结果为" +
                                "\n\t前起落架重量：" + str(result[0]) + "kg" +
                                "\n\t主起落架重量：" + str(result[1]) + "kg" +
                                "\n\t起落架重量：" + str(result[2]) + "kg" +
                                "\n")

            self.history.append([0, input1, input2, result[0], result[1], result[2]])
            print(len(self.history))
            self.nowHistory = len(self.history) - 1

        elif self.choice == 23:
            try:
                input1 = float(self.inputText[0].get())
            except:
                self.logText.insert(1.0, "错误：第一个参数存在格式问题，请检查\n")
                return
            try:
                input2 = float(self.inputText[1].get())
            except:
                self.logText.insert(1.0, "错误：第二个参数存在格式问题，请检查\n")
                return
            try:
                input3 = float(self.inputText[2].get())
            except:
                self.logText.insert(1.0, "错误：第三个参数存在格式问题，请检查\n")
                return

            result = get_result_PW([input1, input2, input3])

            self.outputText[0].delete(1.0, END)
            self.outputText[0].insert(1.0, result[0])
            self.outputText[1].delete(1.0, END)
            self.outputText[1].insert(1.0, result[1])
            self.outputText[2].delete(1.0, END)
            self.outputText[2].insert(1.0, result[2])
            self.logText.insert(1.0, "----------成功----------\n项目：重量分析-尾翼重量估算\n参数为" +
                                "\n\t平尾面积: " + str(input1) + "m²" +
                                "\n\t垂尾面积: " + str(input2) + "m²" +
                                "\n\t设计巡航速度: " + str(input3) + "km/h" +
                                "\n结果为" +
                                "\n\t平尾重量：" + str(result[0]) + "kg" +
                                "\n\t垂尾重量：" + str(result[1]) + "kg" +
                                "\n\t尾翼重量：" + str(result[2]) + "kg" +
                                "\n")
            self.history.append([0, input1, input2, input3, result[0], result[1], result[2]])
            print(len(self.history))
            self.nowHistory = len(self.history) - 1

        elif self.choice == 25:
            try:
                input1 = float(self.inputText[0].get())
            except:
                self.logText.insert(1.0, "错误：第一个参数存在格式问题，请检查\n")
                return
            try:
                input2 = float(self.inputText[1].get())
            except:
                self.logText.insert(1.0, "错误：第二个参数存在格式问题，请检查\n")
                return
            try:
                input3 = float(self.inputText[2].get())
            except:
                self.logText.insert(1.0, "错误：第三个参数存在格式问题，请检查\n")
                return

            result = get_result_FCLQ([input1, input2,input3])

            self.outputText[0].delete(1.0, END)
            self.outputText[0].insert(1.0, result[0])
            self.logText.insert(1.0, "----------成功----------\n项目：结构设计-过载系数-俯冲拉起\n参数为" +
                                "\n\t飞行速度: " + str(input1) + "m/s" +
                                "\n\t重力加速度: " + str(input2) + "m/s²" +
                                "\n\t俯冲半径: " + str(input3) + "m" +
                                "\n结果为" +
                                "\n\t过载系数：" + str(result[0]) + "" +
                                "\n")

            self.history.append([0, input1, input2, input3, result[0]])
            print(len(self.history))
            self.nowHistory = len(self.history) - 1

        elif self.choice == 26:
            try:
                input1 = float(self.inputText[0].get())
            except:
                self.logText.insert(1.0, "错误：第一个参数存在格式问题，请检查\n")
                return

            result = get_result_DSPX([input1])

            self.outputText[0].delete(1.0, END)
            self.outputText[0].insert(1.0, result[0])
            self.logText.insert(1.0, "----------成功----------\n项目：结构设计-过载系数-等速盘旋\n参数为" +
                                "\n\t盘旋倾斜角: " + str(input1) + "°" +
                                "\n结果为" +
                                "\n\t过载系数：" + str(result[0]) + "" +
                                "\n")

            self.history.append([0, input1, result[0]])
            print(len(self.history))
            self.nowHistory = len(self.history) - 1

        elif self.choice == 27:
            try:
                input1 = float(self.inputText[0].get())
            except:
                self.logText.insert(1.0, "错误：第一个参数存在格式问题，请检查\n")
                return
            try:
                input2 = float(self.inputText[1].get())
            except:
                self.logText.insert(1.0, "错误：第二个参数存在格式问题，请检查\n")
                return
            try:
                input3 = float(self.inputText[2].get())
            except:
                self.logText.insert(1.0, "错误：第三个参数存在格式问题，请检查\n")
                return
            try:
                input4 = float(self.inputText[3].get())
            except:
                self.logText.insert(1.0, "错误：第四个参数存在格式问题，请检查\n")
                return
            try:
                input5 = float(self.inputText[4].get())
            except:
                self.logText.insert(1.0, "错误：第五个参数存在格式问题，请检查\n")
                return
            try:
                input6 = float(self.inputText[5].get())
            except:
                self.logText.insert(1.0, "错误：第六个参数存在格式问题，请检查\n")
                return
            try:
                input7 = float(self.inputText[6].get())
            except:
                self.logText.insert(1.0, "错误：第七个参数存在格式问题，请检查\n")
                return

            result = get_result_CZTF([input1, input2, input3, input4, input5, input6, input7])

            self.outputText[0].delete(1.0, END)
            self.outputText[0].insert(1.0, result[0])
            self.outputText[1].delete(1.0, END)
            self.outputText[1].insert(1.0, result[1])
            self.logText.insert(1.0, "----------成功----------\n项目：结构设计-过载系数-垂直突风\n参数为" +
                                "\n\t翼载荷: " + str(input1) + "N/m²" +
                                "\n\t平均气动弦长: " + str(input2) + "m" +
                                "\n\t升力线斜率: " + str(input3) + "/rad" +
                                "\n\t重力加速度: " + str(input4) + "m/s²" +
                                "\n\t飞行高度上的空气密度: " + str(input5) + "kg/m³" +
                                "\n\t飞行速度: " + str(input6) + "m/s" +
                                "\n\t垂直突风速度: " + str(input7) + "m/s" +
                                "\n结果为" +
                                "\n\t离散阵风减缓因子：" + str(result[0]) + "" +
                                "\n\t过载系数：" + str(result[1]) + "" +
                                "\n")

            self.history.append([0, input1, input2, input3, input4,input5,input6,input7, result[0],result[1]])
            print(len(self.history))
            self.nowHistory = len(self.history) - 1

        elif self.choice == 28:
            try:
                input1 = float(self.inputText[0].get())
            except:
                self.logText.insert(1.0, "错误：第一个参数存在格式问题，请检查\n")
                return
            try:
                input2 = float(self.inputText[1].get())
            except:
                self.logText.insert(1.0, "错误：第二个参数存在格式问题，请检查\n")
                return
            try:
                input3 = float(self.inputText[2].get())
            except:
                self.logText.insert(1.0, "错误：第三个参数存在格式问题，请检查\n")
                return
            try:
                input4 = float(self.inputText[3].get())
            except:
                self.logText.insert(1.0, "错误：第四个参数存在格式问题，请检查\n")
                return
            try:
                input5 = float(self.inputText[4].get())
            except:
                self.logText.insert(1.0, "错误：第五个参数存在格式问题，请检查\n")
                return
            try:
                input6 = float(self.inputText[5].get())
            except:
                self.logText.insert(1.0, "错误：第六个参数存在格式问题，请检查\n")
                return

            result = get_result_WDX([input1,input2,input3,input4,input5,input6])

            self.outputText[0].delete(1.0, END)
            self.outputText[0].insert(1.0, result[0])
            self.outputText[1].delete(1.0, END)
            self.outputText[1].insert(1.0, result[1])
            self.logText.insert(1.0, "----------成功----------\n项目：结构设计-结构校核-翼梁稳定性校核\n参数为" +
                                "\n\t剪力: " + str(input1) + "N" +
                                "\n\t上下缘条形心间距离: " + str(input2) + "m" +
                                "\n\t腹板厚度: " + str(input3) + "m" +
                                "\n\t支柱间距: " + str(input4) + "m" +
                                "\n\t弹性模量: " + str(input5) + "N/m²" +
                                "\n\t泊松比: " + str(input6) + "" +
                                "\n结果为" +
                                "\n\t名义剪应力：" + str(result[0]) + "N/m²" +
                                "\n\t临界失稳剪应力：" + str(result[1]) + "N/m²" +
                                "\n")

            self.history.append([0, input1,input2,input3,input4,input5,input6, result[0],result[1]])
            print(len(self.history))
            self.nowHistory = len(self.history) - 1

        elif self.choice == 29:
            try:
                input1 = float(self.inputText[0].get())
            except:
                self.logText.insert(1.0, "错误：第一个参数存在格式问题，请检查\n")
                return
            try:
                input2 = float(self.inputText[1].get())
            except:
                self.logText.insert(1.0, "错误：第二个参数存在格式问题，请检查\n")
                return
            try:
                input3 = float(self.inputText[2].get())
            except:
                self.logText.insert(1.0, "错误：第三个参数存在格式问题，请检查\n")
                return

            result = get_result_QD([input1, input2,input3])

            self.outputText[0].delete(1.0, END)
            self.outputText[0].insert(1.0, result[0])
            self.logText.insert(1.0, "----------成功----------\n项目：结构设计-结构校核-翼梁强度校核\n参数为" +
                                "\n\t弯矩: " + str(input1) + "N*m" +
                                "\n\t截面最大高度: " + str(input2) + "m" +
                                "\n\t截面惯性矩: " + str(input3) + "kg*m²" +
                                "\n结果为" +
                                "\n\t最大正应力：" + str(result[0]) + "N/m²" +
                                "\n")
            self.history.append([0, input1, input2,input3, result[0]])
            print(len(self.history))
            self.nowHistory = len(self.history) - 1

        elif self.choice == 30:
            try:
                input1 = float(self.inputText[0].get())
            except:
                self.logText.insert(1.0, "错误：第一个参数存在格式问题，请检查\n")
                return
            try:
                input2 = float(self.inputText[1].get())
            except:
                self.logText.insert(1.0, "错误：第二个参数存在格式问题，请检查\n")
                return
            try:
                input3 = float(self.inputText[2].get())
            except:
                self.logText.insert(1.0, "错误：第三个参数存在格式问题，请检查\n")
                return
            try:
                input4 = float(self.inputText[3].get())
            except:
                self.logText.insert(1.0, "错误：第四个参数存在格式问题，请检查\n")
                return

            result = get_result_GZB_TZB([input1,input2,input3,input4])

            self.outputText[0].delete(1.0, END)
            self.outputText[0].insert(1.0, result[0])
            self.logText.insert(1.0, "----------成功----------\n项目：总体设计-功重比→推重比\n参数为" +
                                "\n\t功重比: " + str(input1) + "hp/kg" +
                                "\n\t速度: " + str(input2) + "m/s" +
                                "\n\t发动机效率: " + str(input3) + " " +
                                "\n\t高度修正因子: " + str(input4) + " " +
                                "\n结果为" +
                                "\n\t推重比：" + str(result[0]) + "N/kg" +
                                "\n")

            self.history.append([0, input1,input2,input3,input4,result[0]])
            print(len(self.history))
            self.nowHistory = len(self.history) - 1

        elif self.choice == 31:
            try:
                input1 = float(self.inputText[0].get())
            except:
                self.logText.insert(1.0, "错误：第一个参数存在格式问题，请检查\n")
                return
            try:
                input2 = float(self.inputText[1].get())
            except:
                self.logText.insert(1.0, "错误：第二个参数存在格式问题，请检查\n")
                return
            try:
                input3 = float(self.inputText[2].get())
            except:
                self.logText.insert(1.0, "错误：第三个参数存在格式问题，请检查\n")
                return
            try:
                input4 = float(self.inputText[3].get())
            except:
                self.logText.insert(1.0, "错误：第四个参数存在格式问题，请检查\n")
                return

            result = get_result_Gas_Ferrar([input1,input2,input3,input4])

            self.outputText[0].delete(1.0, END)
            self.outputText[0].insert(1.0, result[0])
            self.logText.insert(1.0, "----------成功----------\n项目：总体设计-推重比→功重比\n参数为" +
                                "\n\t推重比: " + str(input1) + "N/kg" +
                                "\n\t速度: " + str(input2) + "m/s" +
                                "\n\t发动机效率: " + str(input3) + " " +
                                "\n\t高度修正因子: " + str(input4) + " " +
                                "\n结果为" +
                                "\n\t功重比：" + str(result[0]) + "hp/kg" +
                                "\n")

            self.history.append([0, input1,input2,input3,input4,result[0]])
            print(len(self.history))
            self.nowHistory = len(self.history) - 1

        elif self.choice == 32:
            try:
                input1 = float(self.inputText[0].get())
            except:
                self.logText.insert(1.0, "错误：第一个参数存在格式问题，请检查\n")
                return
            try:
                input2 = float(self.inputText[1].get())
            except:
                self.logText.insert(1.0, "错误：第二个参数存在格式问题，请检查\n")
                return
            try:
                input3 = float(self.inputText[2].get())
            except:
                self.logText.insert(1.0, "错误：第三个参数存在格式问题，请检查\n")
                return
            try:
                input4 = float(self.inputText[3].get())
            except:
                self.logText.insert(1.0, "错误：第四个参数存在格式问题，请检查\n")
                return

            result = get_result_HLJ([input1, input2, input3, input4])

            self.outputText[0].delete(1.0, END)
            self.outputText[0].insert(1.0, result[0])
            self.logText.insert(1.0, "----------成功----------\n项目：气动设计-机翼-1/n弦线的后掠角\n参数为" +
                                "\n\t前缘后掠角: " + str(input1) + "°" +
                                "\n\tn: " + str(input2) + "" +
                                "\n\t展弦比: " + str(input3) + " " +
                                "\n\t根梢比: " + str(input4) + " " +
                                "\n结果为" +
                                "\n\t1/n弦线的后掠角：" + str(result[0]) + "°" +
                                "\n")

            self.history.append([0, input1, input2, input3, input4, result[0]])
            print(len(self.history))
            self.nowHistory = len(self.history) - 1

        elif self.choice == 33:
            try:
                input1 = float(self.inputText[0].get())
            except:
                self.logText.insert(1.0, "错误：第一个参数存在格式问题，请检查\n")
                return
            try:
                input2 = float(self.inputText[1].get())
            except:
                self.logText.insert(1.0, "错误：第二个参数存在格式问题，请检查\n")
                return
            try:
                input3 = float(self.inputText[2].get())
            except:
                self.logText.insert(1.0, "错误：第三个参数存在格式问题，请检查\n")
                return
            try:
                input4 = float(self.inputText[3].get())
            except:
                self.logText.insert(1.0, "错误：第四个参数存在格式问题，请检查\n")
                return

            result = get_result_SLXXL([input1, input2, input3, input4])

            self.outputText[0].delete(1.0, END)
            self.outputText[0].insert(1.0, result[0])
            self.outputText[1].delete(1.0, END)
            self.outputText[1].insert(1.0, result[1])
            self.logText.insert(1.0, "----------成功----------\n项目：性能分析-升阻性能-机翼升力线斜率\n参数为" +
                                "\n\t机翼展弦比: " + str(input1) + "" +
                                "\n\t机翼前缘后掠角: " + str(input2) + "°" +
                                "\n\t机翼后缘后掠角: " + str(input3) + "°" +
                                "\n\t机翼根梢比: " + str(input4) + " " +
                                "\n结果为" +
                                "\n\t机翼升力线斜率：" + str(result[0]) + "/rad" +
                                "\n\t机翼升力线斜率(*pi)：" + str(result[1]) + "/rad" 
                                "\n")

            self.history.append([0, input1, input2, input3, input4, result[0],result[1]])
            print(len(self.history))
            self.nowHistory = len(self.history) - 1

        elif self.choice == 34:
            try:
                input1 = float(self.inputText[0].get())
            except:
                self.logText.insert(1.0, "错误：第一个参数存在格式问题，请检查\n")
                return

            result = get_result_Oswald([input1])

            self.outputText[0].delete(1.0, END)
            self.outputText[0].insert(1.0, result[0])
            self.logText.insert(1.0, "----------成功----------\n项目：性能分析-升阻性能-奥斯瓦尔德效率因子\n参数为" +
                                "\n\t机翼展弦比: " + str(input1) + "" +
                                "\n结果为" +
                                "\n\t奥斯瓦尔德效率因子：" + str(result[0]) + "" +
                                "\n")

            self.history.append([0, input1, result[0]])
            print(len(self.history))
            self.nowHistory = len(self.history) - 1

        elif self.choice == 35:
            try:
                input1 = float(self.inputText[0].get())
            except:
                self.logText.insert(1.0, "错误：第一个参数存在格式问题，请检查\n")
                return
            try:
                input2 = float(self.inputText[1].get())
            except:
                self.logText.insert(1.0, "错误：第二个参数存在格式问题，请检查\n")
                return
            try:
                input3 = float(self.inputText[2].get())
            except:
                self.logText.insert(1.0, "错误：第三个参数存在格式问题，请检查\n")
                return

            result = get_result_PF([input1,input2,input3])

            self.outputText[0].delete(1.0, END)
            self.outputText[0].insert(1.0, result[0])
            self.logText.insert(1.0, "----------成功----------\n项目：性能分析-升阻性能-巡航升力系数\n参数为" +
                                "\n\t翼载荷: " + str(input1) + "kg/m²" +
                                "\n\t当地空气密度: " + str(input2) + "kg/m³" +
                                "\n\t巡航速度: " + str(input3) + "m/s" +
                                "\n结果为" +
                                "\n\t巡航升力系数：" + str(result[0]) + "" +
                                "\n")

            self.history.append([0, input1,input2,input3, result[0]])
            print(len(self.history))
            self.nowHistory = len(self.history) - 1

        elif self.choice == 36:
            try:
                input1 = float(self.inputText[0].get())
            except:
                self.logText.insert(1.0, "错误：第一个参数存在格式问题，请检查\n")
                return
            try:
                input2 = float(self.inputText[1].get())
            except:
                self.logText.insert(1.0, "错误：第二个参数存在格式问题，请检查\n")
                return
            try:
                input3 = float(self.inputText[2].get())
            except:
                self.logText.insert(1.0, "错误：第三个参数存在格式问题，请检查\n")
                return

            result = get_result_DI([input1, input2, input3])

            self.outputText[0].delete(1.0, END)
            self.outputText[0].insert(1.0, result[0])
            self.logText.insert(1.0, "----------成功----------\n项目：性能分析-升阻性能-诱导阻力系数\n参数为" +
                                "\n\t升力系数: " + str(input1) + "" +
                                "\n\t奥斯瓦尔德效率因子: " + str(input2) + "" +
                                "\n\t机翼展弦比: " + str(input3) + "" +
                                "\n结果为" +
                                "\n\t诱导阻力系数：" + str(result[0]) + "" +
                                "\n")

            self.history.append([0, input1, input2, input3, result[0]])
            print(len(self.history))
            self.nowHistory = len(self.history) - 1

        elif self.choice == 37:
            try:
                input1 = float(self.inputText[0].get())
            except:
                self.logText.insert(1.0, "错误：第一个参数存在格式问题，请检查\n")
                return
            try:
                input2 = float(self.inputText[1].get())
            except:
                self.logText.insert(1.0, "错误：第二个参数存在格式问题，请检查\n")
                return
            try:
                input3 = float(self.inputText[2].get())
            except:
                self.logText.insert(1.0, "错误：第三个参数存在格式问题，请检查\n")
                return
            try:
                input4 = float(self.inputText[3].get())
            except:
                self.logText.insert(1.0, "错误：第四个参数存在格式问题，请检查\n")
                return
            try:
                input5 = float(self.inputText[4].get())
            except:
                self.logText.insert(1.0, "错误：第五个参数存在格式问题，请检查\n")
                return
            try:
                input6 = float(self.inputText[5].get())
            except:
                self.logText.insert(1.0, "错误：第六个参数存在格式问题，请检查\n")
                return
            try:
                input7 = float(self.inputText[6].get())
            except:
                self.logText.insert(1.0, "错误：第七个参数存在格式问题，请检查\n")
                return
            try:
                input8 = float(self.inputText[7].get())
            except:
                self.logText.insert(1.0, "错误：第八个参数存在格式问题，请检查\n")
                return
            try:
                input9 = float(self.inputText[8].get())
            except:
                self.logText.insert(1.0, "错误：第九个参数存在格式问题，请检查\n")
                return
            try:
                input10 = float(self.inputText[9].get())
            except:
                self.logText.insert(1.0, "错误：第十个参数存在格式问题，请检查\n")
                return

            result = get_result_TAKEOFF_distance([input1, input2, input3, input4, input5, input6, input7, input8, input9, input10])

            self.outputText[0].delete(1.0, END)
            self.outputText[0].insert(1.0, result[0])
            self.logText.insert(1.0, "----------成功----------\n项目：性能分析-起降性能-起飞滑跑距离\n参数为" +
                                "\n\t起飞重量(N): " + str(input1) + "N" +
                                "\n\t重力加速度: " + str(input2) + "m/s²" +
                                "\n\t机翼面积: " + str(input3) + "m²" +
                                "\n\t当地空气密度: " + str(input4) + "kg/m³" +
                                "\n\t最大升力系数: " + str(input5) + "" +
                                "\n\t发动机推力: " + str(input6) + "N" +
                                "\n\t地面摩擦阻力系数: " + str(input7) + "" +
                                "\n\t零升阻力系数: " + str(input8) + "" +
                                "\n\t奥斯瓦尔德效率因子: " + str(input9) + "" +
                                "\n\t机翼展弦比: " + str(input10) + "" +
                                "\n结果为" +
                                "\n\t起飞滑跑距离：" + str(result[0]) + "m" +
                                "\n")

            self.history.append([0, input1, input2, input3, input4, input5, input6, input7, input8, input9, input10, result[0]])
            print(len(self.history))
            self.nowHistory = len(self.history) - 1

        elif self.choice == 38:
            try:
                input1 = float(self.inputText[0].get())
            except:
                self.logText.insert(1.0, "错误：第一个参数存在格式问题，请检查\n")
                return
            try:
                input2 = float(self.inputText[1].get())
            except:
                self.logText.insert(1.0, "错误：第二个参数存在格式问题，请检查\n")
                return
            try:
                input3 = float(self.inputText[2].get())
            except:
                self.logText.insert(1.0, "错误：第三个参数存在格式问题，请检查\n")
                return

            result = get_result_TAKEOFF_velocity([input1, input2, input3])

            self.outputText[0].delete(1.0, END)
            self.outputText[0].insert(1.0, result[0])
            self.logText.insert(1.0, "----------成功----------\n项目：性能分析-起降性能-起飞离地速度\n参数为" +
                                "\n\t翼载荷: " + str(input1) + "kg/m²" +
                                "\n\t当地空气密度: " + str(input2) + "kg/m³" +
                                "\n\t最大升力系数: " + str(input3) + "" +
                                "\n结果为" +
                                "\n\t起飞离地速度：" + str(result[0]) + "m/s" +
                                "\n")

            self.history.append([0, input1, input2, input3, result[0]])
            print(len(self.history))
            self.nowHistory = len(self.history) - 1

        elif self.choice == 39:
            try:
                input1 = float(self.inputText[0].get())
            except:
                self.logText.insert(1.0, "错误：第一个参数存在格式问题，请检查\n")
                return
            try:
                input2 = float(self.inputText[1].get())
            except:
                self.logText.insert(1.0, "错误：第二个参数存在格式问题，请检查\n")
                return
            try:
                input3 = float(self.inputText[2].get())
            except:
                self.logText.insert(1.0, "错误：第三个参数存在格式问题，请检查\n")
                return

            result = get_result_LAND_velocity([input1, input2, input3])

            self.outputText[0].delete(1.0, END)
            self.outputText[0].insert(1.0, result[0])
            self.logText.insert(1.0, "----------成功----------\n项目：性能分析-起降性能-进场速度\n参数为" +
                                "\n\t翼载荷: " + str(input1) + "kg/m²" +
                                "\n\t当地空气密度: " + str(input2) + "kg/m³" +
                                "\n\t最大升力系数: " + str(input3) + "" +
                                "\n结果为" +
                                "\n\t起飞离地速度：" + str(result[0]) + "m/s" +
                                "\n")

            self.history.append([0, input1, input2, input3, result[0]])
            print(len(self.history))
            self.nowHistory = len(self.history) - 1

        elif self.choice == 40:
            try:
                input1 = float(self.inputText[0].get())
            except:
                self.logText.insert(1.0, "错误：第一个参数存在格式问题，请检查\n")
                return
            try:
                input2 = float(self.inputText[1].get())
            except:
                self.logText.insert(1.0, "错误：第二个参数存在格式问题，请检查\n")
                return
            try:
                input3 = float(self.inputText[2].get())
            except:
                self.logText.insert(1.0, "错误：第三个参数存在格式问题，请检查\n")
                return
            try:
                input4 = float(self.inputText[3].get())
            except:
                self.logText.insert(1.0, "错误：第四个参数存在格式问题，请检查\n")
                return
            try:
                input5 = float(self.inputText[4].get())
            except:
                self.logText.insert(1.0, "错误：第五个参数存在格式问题，请检查\n")
                return

            result = get_result_Cruise_distance([input1, input2, input3, input4, input5])

            self.outputText[0].delete(1.0, END)
            self.outputText[0].insert(1.0, result[0])
            self.logText.insert(1.0, "----------成功----------\n项目：性能分析-续航性能-等速巡航航程\n参数为" +
                                "\n\t巡航速度: " + str(input1) + "m/s" +
                                "\n\t耗油率: " + str(input2) + "N/(kW·h)" +
                                "\n\t升阻比: " + str(input3) + "" +
                                "\n\t巡航段初始重量: " + str(input4) + "N" +
                                "\n\t巡航段结束重量: " + str(input5) + "N" 
                                "\n结果为" +
                                "\n\t等速巡航航程：" + str(result[0]) + "m" +
                                "\n")

            self.history.append([0, input1, input2, input3, input4, input5, result[0]])
            print(len(self.history))
            self.nowHistory = len(self.history) - 1

        elif self.choice == 41:
            try:
                input1 = float(self.inputText[0].get())
            except:
                self.logText.insert(1.0, "错误：第一个参数存在格式问题，请检查\n")
                return
            try:
                input2 = float(self.inputText[1].get())
            except:
                self.logText.insert(1.0, "错误：第二个参数存在格式问题，请检查\n")
                return
            try:
                input3 = float(self.inputText[2].get())
            except:
                self.logText.insert(1.0, "错误：第三个参数存在格式问题，请检查\n")
                return
            try:
                input4 = float(self.inputText[3].get())
            except:
                self.logText.insert(1.0, "错误：第四个参数存在格式问题，请检查\n")
                return

            result = get_result_Cruise_time([input1, input2, input3, input4])

            self.outputText[0].delete(1.0, END)
            self.outputText[0].insert(1.0, result[0])
            self.logText.insert(1.0, "----------成功----------\n项目：性能分析-续航性能-等速巡航航程\n参数为" +
                                "\n\t耗油率: " + str(input1) + "N/(kW·h)" +
                                "\n\t升阻比: " + str(input2) + "" +
                                "\n\t巡航段初始重量: " + str(input3) + "N" +
                                "\n\t巡航段结束重量: " + str(input4) + "N" +
                                "\n结果为" +
                                "\n\t航时：" + str(result[0]) + "h" +
                                "\n")

            self.history.append([0, input1, input2, input3, input4, result[0]])
            print(len(self.history))
            self.nowHistory = len(self.history) - 1

        elif self.choice == 42:
            try:
                input1 = float(self.inputText[0].get())
            except:
                self.logText.insert(1.0, "错误：第一个参数存在格式问题，请检查\n")
                return
            try:
                input2 = float(self.inputText[1].get())
            except:
                self.logText.insert(1.0, "错误：第二个参数存在格式问题，请检查\n")
                return
            try:
                input3 = float(self.inputText[2].get())
            except:
                self.logText.insert(1.0, "错误：第三个参数存在格式问题，请检查\n")
                return

            result = get_result_hover([input1, input2, input3])

            self.outputText[0].delete(1.0, END)
            self.outputText[0].insert(1.0, result[0])
            self.outputText[1].delete(1.0, END)
            self.outputText[1].insert(1.0, result[1])
            self.outputText[2].delete(1.0, END)
            self.outputText[2].insert(1.0, result[2])
            self.outputText[3].delete(1.0, END)
            self.outputText[3].insert(1.0, result[3])
            self.logText.insert(1.0, "----------成功----------\n项目：性能分析-盘旋性能-盘旋过载&盘旋角速率&盘旋半径&盘旋一周所需时间\n参数为" +
                                "\n\t重力加速度: " + str(input1) + "m/s²" +
                                "\n\t盘旋过载: " + str(input2) + "" +
                                "\n\t飞行速度: " + str(input3) + "m/s" +
                                "\n结果为" +
                                "\n\t盘旋盘旋过载：" + str(result[0]) + "" +
                                "\n\t盘旋角速率：" + str(result[1]) + "rad/s" +
                                "\n\t盘旋半径：" + str(result[2]) + "m" +
                                "\n\t盘旋一周所需时间：" + str(result[3]) + "s" +
                                "\n")

            self.history.append([0, input1, input2, input3,result[0],result[1],result[2]])
            print(len(self.history))
            self.nowHistory = len(self.history) - 1

        elif self.choice == 43:
            try:
                input1 = float(self.inputText[0].get())
            except:
                self.logText.insert(1.0, "错误：第一个参数存在格式问题，请检查\n")
                return
            try:
                input2 = float(self.inputText[1].get())
            except:
                self.logText.insert(1.0, "错误：第二个参数存在格式问题，请检查\n")
                return
            try:
                input3 = float(self.inputText[2].get())
            except:
                self.logText.insert(1.0, "错误：第三个参数存在格式问题，请检查\n")
                return

            result = get_result_E_max([input1, input2, input3])

            self.outputText[0].delete(1.0, END)
            self.outputText[0].insert(1.0, result[0])
            self.logText.insert(1.0, "----------成功----------\n项目：性能分析-升阻性能-最大升阻比\n参数为" +
                                "\n\t零升阻力系数: " + str(input1) + "" +
                                "\n\t奥斯瓦尔德效率因子: " + str(input2) + "" +
                                "\n\t展弦比: " + str(input3) + "" +
                                "\n结果为" +
                                "\n\t最大升阻比：" + str(result[0]) + "" +
                                "\n")

            self.history.append([0, input1, input2, input3, result[0]])
            print(len(self.history))
            self.nowHistory = len(self.history) - 1

        elif self.choice == 44:
            try:
                input1 = float(self.inputText[0].get())
            except:
                self.logText.insert(1.0, "错误：第一个参数存在格式问题，请检查\n")
                return
            try:
                input2 = float(self.inputText[1].get())
            except:
                self.logText.insert(1.0, "错误：第二个参数存在格式问题，请检查\n")
                return
            try:
                input3 = float(self.inputText[2].get())
            except:
                self.logText.insert(1.0, "错误：第三个参数存在格式问题，请检查\n")
                return

            result = get_result_H_n([input1, input2, input3])

            self.outputText[0].delete(1.0, END)
            self.outputText[0].insert(1.0, result[0])
            self.logText.insert(1.0, "----------成功----------\n项目：稳定性分析-纵向静稳定裕度\n参数为" +
                                "\n\t全机焦点位置: " + str(input1) + "m" +
                                "\n\t重心位置: " + str(input2) + "m" +
                                "\n\t平均气动弦长: " + str(input3) + "m" +
                                "\n结果为" +
                                "\n\t纵向静稳定裕度：" + str(result[0]) + "" +
                                "\n")

            self.history.append([0, input1, input2, input3, result[0]])
            print(len(self.history))
            self.nowHistory = len(self.history) - 1

        elif self.choice == 45:
            try:
                input1 = float(self.inputText[0].get())
            except:
                self.logText.insert(1.0, "错误：第一个参数存在格式问题，请检查\n")
                return
            try:
                input2 = float(self.inputText[1].get())
            except:
                self.logText.insert(1.0, "错误：第二个参数存在格式问题，请检查\n")
                return
            try:
                input3 = float(self.inputText[2].get())
            except:
                self.logText.insert(1.0, "错误：第三个参数存在格式问题，请检查\n")
                return

            result = get_result_air([input1, input2, input3])

            self.outputText[0].delete(1.0, END)
            self.outputText[0].insert(1.0, result[0])
            self.outputText[1].delete(1.0, END)
            self.outputText[1].insert(1.0, result[1])
            self.outputText[2].delete(1.0, END)
            self.outputText[2].insert(1.0, result[2])
            self.outputText[3].delete(1.0, END)
            self.outputText[3].insert(1.0, result[3])
            self.outputText[4].delete(1.0, END)
            self.outputText[4].insert(1.0, result[4])
            self.outputText[5].delete(1.0, END)
            self.outputText[5].insert(1.0, result[5])
            self.outputText[6].delete(1.0, END)
            self.outputText[6].insert(1.0, result[6])
            self.outputText[7].delete(1.0, END)
            self.outputText[7].insert(1.0, result[7])
            self.outputText[8].delete(1.0, END)
            self.outputText[8].insert(1.0, result[8])
            self.outputText[9].delete(1.0, END)
            self.outputText[9].insert(1.0, result[9])
            self.logText.insert(1.0, "----------成功----------\n项目：空气动力学-大气参数\n参数为" +
                                "\n\t高度: " + str(input1) + "m" +
                                "\n\t速度: " + str(input2) + "m/s" +
                                "\n\t特征尺寸: " + str(input3) + "m" +
                                "\n结果为" +
                                "\n\t温度(K)：" + str(result[0]) + "K" +
                                "\n\t温度(摄氏度)：" + str(result[1]) + "℃" +
                                "\n\t密度：" + str(result[2]) + "kg/m³" +
                                "\n\t气压：" + str(result[3]) + "Pa" +
                                "\n\t动压：" + str(result[4]) + "N/m²" +
                                "\n\t音速：" + str(result[5]) + "m/s" +
                                "\n\t粘度：" + str(result[6]) + "N*s/m²" +
                                "\n\t马赫数：" + str(result[7]) + "" +
                                "\n\t雷诺数：" + str(result[8]) + "" +
                                "\n\t密度比：" + str(result[9]) + " " +
                                "\n")

            self.history.append([0, input1, input2, input3, result[0], result[1], result[2], result[3], result[4], result[5], result[6], result[7], result[8], result[9]])
            print(len(self.history))
            self.nowHistory = len(self.history) - 1

        elif self.choice == 46:
            try:
                input1 = float(self.inputText[0].get())
            except:
                self.logText.insert(1.0, "错误：第一个参数存在格式问题，请检查\n")
                return
            try:
                input2 = float(self.inputText[1].get())
            except:
                self.logText.insert(1.0, "错误：第二个参数存在格式问题，请检查\n")
                return
            try:
                input3 = float(self.inputText[2].get())
            except:
                self.logText.insert(1.0, "错误：第三个参数存在格式问题，请检查\n")
                return

            result = get_result_lift_airfoil([input1, input2,input3])

            self.outputText[0].delete(1.0, END)
            self.outputText[0].insert(1.0, result[0])
            self.logText.insert(1.0, "----------成功----------\n项目：空气动力学-翼型升力\n参数为" +
                                "\n\t空气密度: " + str(input1) + "kg/m³" +
                                "\n\t来流速度: " + str(input2) + "m/s" +
                                "\n\t速度环量: " + str(input3) + "m²/s" +
                                "\n结果为" +
                                "\n\t翼型升力：" + str(result[0]) + "N/m" +
                                "\n")
            self.history.append([0, input1, input2,input3, result[0]])
            print(len(self.history))
            self.nowHistory = len(self.history) - 1

        elif self.choice == 47:
            try:
                input1 = float(self.inputText[0].get())
            except:
                self.logText.insert(1.0, "错误：第一个参数存在格式问题，请检查\n")
                return
            try:
                input2 = float(self.inputText[1].get())
            except:
                self.logText.insert(1.0, "错误：第二个参数存在格式问题，请检查\n")
                return
            try:
                input3 = float(self.inputText[2].get())
            except:
                self.logText.insert(1.0, "错误：第三个参数存在格式问题，请检查\n")
                return

            result = get_result_bonuli([input1, input2, input3])

            self.outputText[0].delete(1.0, END)
            self.outputText[0].insert(1.0, result[0])
            self.logText.insert(1.0, "----------成功----------\n项目：空气动力学-伯努利方程(皮托管测空速)\n参数为" +
                                "\n\t总压(皮托管测得压强): " + str(input1) + "Pa" +
                                "\n\t大气静压: " + str(input2) + "Pa" +
                                "\n\t大气密度: " + str(input3) + "kg/m³" +
                                "\n结果为" +
                                "\n\t飞行速度：" + str(result[0]) + "m/s" +
                                "\n")
            self.history.append([0, input1, input2, input3, result[0]])
            print(len(self.history))
            self.nowHistory = len(self.history) - 1

        elif self.choice == 48:
            try:
                input1 = float(self.inputText[0].get())
            except:
                self.logText.insert(1.0, "错误：第一个参数存在格式问题，请检查\n")
                return

            result = get_result_mahe_angle([input1])

            self.outputText[0].delete(1.0, END)
            self.outputText[0].insert(1.0, result[0])
            self.logText.insert(1.0, "----------成功----------\n项目：空气动力学-马赫角\n参数为" +
                                "\n\t马赫数: " + str(input1) + "" +
                                "\n结果为" +
                                "\n\t马赫角：" + str(result[0]) + "°" +
                                "\n")
            self.history.append([0, input1, result[0]])
            print(len(self.history))
            self.nowHistory = len(self.history) - 1

        elif self.choice == 49:
            try:
                input1 = float(self.inputText[0].get())
            except:
                self.logText.insert(1.0, "错误：第一个参数存在格式问题，请检查\n")
                return
            try:
                input2 = float(self.inputText[1].get())
            except:
                self.logText.insert(1.0, "错误：第二个参数存在格式问题，请检查\n")
                return

            result = get_result_linear([input1, input2])

            self.outputText[0].delete(1.0, END)
            self.outputText[0].insert(1.0, result[0])
            self.outputText[1].delete(1.0, END)
            self.outputText[1].insert(1.0, result[1])
            self.logText.insert(1.0, "----------成功----------\n项目：空气动力学-线化超音速理论(小迎角平板翼型)\n参数为" +
                                "\n\t迎角: " + str(input1) + "°" +
                                "\n\t马赫数: " + str(input2) + "" +
                                "\n结果为" +
                                "\n\t升力系数：" + str(result[0]) + "" +
                                "\n\t波阻系数：" + str(result[1]) + "" +
                                "\n")
            self.history.append([0, input1, input2, result[0], result[1]])
            print(len(self.history))
            self.nowHistory = len(self.history) - 1

        elif self.choice == 50:
            try:
                input1 = float(self.inputText[0].get())
            except:
                self.logText.insert(1.0, "错误：第一个参数存在格式问题，请检查\n")
                return
            try:
                input2 = float(self.inputText[1].get())
            except:
                self.logText.insert(1.0, "错误：第二个参数存在格式问题，请检查\n")
                return
            try:
                input3 = float(self.inputText[2].get())
            except:
                self.logText.insert(1.0, "错误：第三个参数存在格式问题，请检查\n")
                return
            try:
                input4 = float(self.inputText[3].get())
            except:
                self.logText.insert(1.0, "错误：第四个参数存在格式问题，请检查\n")
                return

            result = get_result_boundary_layer([input1, input2,input3,input4])

            self.outputText[0].delete(1.0, END)
            self.outputText[0].insert(1.0, result[0])
            self.history.append([0, input1, input2, result[0]])
            self.logText.insert(1.0, "----------成功----------\n项目：空气动力学-平板层流边界层厚度\n参数为" +
                                "\n\t距离平板前缘距离: " + str(input1) + "m" +
                                "\n\t气体密度: " + str(input2) + "kg/m³" +
                                "\n\t气体粘度: " + str(input2) + "N*s/m²" +
                                "\n\t来流速度: " + str(input2) + "m/s" +
                                "\n结果为" +
                                "\n\t边界层厚度：" + str(result[0]) + "m" +
                                "\n")
            print(len(self.history))
            self.nowHistory = len(self.history) - 1

    def save(self):
        filename = asksaveasfilename(initialdir='C:/', filetypes=[('FlyDesignFile', '.fdf')])
        print(filename)
        if not str.endswith(filename, ".fdf"):
            try:
                os.rename(filename, filename + '.fdf')
            except:
                print("add index failed")
            filename = filename + '.fdf'
        file = open(filename, 'w')
        for index in range(0, len(self.history)):
            for rowIndex in range(0, len(self.history[index])):
                file.write(str(self.history[index][rowIndex]))
                file.write(' ')
            file.write('\n')

    def read(self):
        filename = askopenfilename(initialdir='C:/', filetypes=[('FlyDesignFile', '.fdf')])
        file = open(filename, 'r')
        now = file.readline()
        self.history = []
        while now:
            numbers = str.split(now, ' ')
            first = int(numbers[0])
            if first == 0:
                input1 = float(numbers[1])
                input2 = float(numbers[2])
                input3 = float(numbers[3])
                input4 = float(numbers[4])
                output1 = float(numbers[5])
                self.history.append([0, input1, input2, input3, input4, output1])
                print(len(self.history))
                self.nowHistory = len(self.history)
            now = file.readline()
        self.last_history()

    def last_history(self):

        if self.nowHistory == -1:
            return
        self.nowHistory -= 1
        if self.nowHistory == -1:
            return
        if self.history[self.nowHistory][0] == 0:
            self.callback_main_qifei()
            for index in range(len(self.inputText)):
                if self.inputText[index] == None:
                    break
                self.inputText[index].insert(0, str(self.history[self.nowHistory][index + 1]))
            for index in range(len(self.outputText)):
                if self.outputText[index] == None:
                    break
                self.outputText[index].insert(1.0, str(self.history[self.nowHistory][index + 5]))

    def next_history(self):
        if self.nowHistory >= len(self.history):
            return
        self.nowHistory += 1
        if self.nowHistory >= len(self.history):
            return
        if self.history[self.nowHistory][0] == 0:
            self.callback_main_qifei()
            for index in range(len(self.inputText)):
                if self.inputText[index] == None:
                    break
                self.inputText[index].insert(0, str(self.history[self.nowHistory][index + 1]))
            for index in range(len(self.outputText)):
                if self.outputText[index] == None:
                    break
                self.outputText[index].insert(1.0, str(self.history[self.nowHistory][index + 5]))

    # 设置窗口
    def set_init_window(self):
        self.history = []
        self.nowHistory = -1
        self.root.title("飞机设计与分析快速计算软件_v1.0")  # 窗口名
        self.root.geometry('1080x720+20+20')  # 窗口位置和偏移量
        self.root["bg"] = "#EEEEEE"  # 窗口背景色
        # self.root.attributes("-alpha",0.9)                          #虚化，值越小虚化程度越高

        menubar = Menu(self.root)  # 顶部菜单
        fmenu = Menu(menubar, tearoff=False)  # 在 menubar菜单实例上建立新的子菜单实例
        fmenu.add_command(label='新建', command=self.save)
        fmenu.add_command(label='打开', command=self.read)
        fmenu.add_command(label='保存', command=self.save)
        fmenu.add_command(label='另存为', command=self.save)

        emenu = Menu(menubar, tearoff=False)  # 在 menubar菜单实例上建立新的子菜单实例
        mainDesign = Menu(emenu, tearoff=False)
        self.choice = 0
        # 总体设计子菜单
        mainDesign.add_command(label='起飞重量', command=self.callback_main_qifei)
        mainDesign.add_command(label='零升阻力系数', command=self.callback_main_linsheng)
        mainDesign.add_command(label='翼载荷', command=self.callback_main_yizai)
        mainDesign.add_command(label='推重比', command=self.callback_main_tuizhong)
        mainDesign.add_command(label='功重比', command=self.callback_main_gongzhong)
        mainDesign.add_command(label='功重比→推重比', command=self.callback_gongzhong2tuizhong)
        mainDesign.add_command(label='推重比→功重比', command=self.callback_tuizhong2gongzhong)

        # 气动设计子菜单
        airDesign = Menu(emenu, tearoff=False)

        wingDesign = Menu(airDesign, tearoff=False)
        wingDesign.add_command(label='机翼展弦比', command=self.callback_air_jiyi_zhanxian)
        wingDesign.add_command(label='机翼根梢比', command=self.callback_air_jiyi_genshao)
        wingDesign.add_command(label='机翼根弦长和梢弦长', command=self.callback_air_jiyi_xianchang)
        wingDesign.add_command(label='机翼平均气动弦长', command=self.callback_air_jiyi_qidong)
        wingDesign.add_command(label='1/n弦线的后掠角', command=self.callback_HLJ)
        airDesign.add_cascade(label="机翼", menu=wingDesign)

        fuselageDesign = Menu(airDesign, tearoff=False)
        fuselageDesign.add_command(label='机身长细比', command=self.callback_air_jishen_changxi)
        fuselageDesign.add_command(label='纵向机身容量参数', command=self.callback_air_pingwei_zongxiang)
        fuselageDesign.add_command(label='航向机身容量参数', command=self.callback_air_chuiwei_canshu)
        airDesign.add_cascade(label="机身", menu=fuselageDesign)

        flatTailDesign = Menu(airDesign, tearoff=False)
        flatTailDesign.add_command(label='平尾尾容量', command=self.callback_air_pingwei_rongliang)
        flatTailDesign.add_command(label='平尾展弦比', command=self.callback_air_pingwei_zhanxian)
        flatTailDesign.add_command(label='平尾根梢比', command=self.callback_air_pingwei_gengshao)
        flatTailDesign.add_command(label='平尾根弦长和梢弦长', command=self.callback_air_pingwei_xianchang)
        airDesign.add_cascade(label="平尾", menu=flatTailDesign)

        verticalTailDesign = Menu(airDesign, tearoff=False)
        verticalTailDesign.add_command(label='垂尾尾容量', command=self.callback_air_chuiwei_rongliang)
        verticalTailDesign.add_command(label='垂尾展弦比', command=self.callback_air_chuiwei_zhanxian)
        verticalTailDesign.add_command(label='垂尾根梢比', command=self.callback_air_chuiwei_gengshao)
        verticalTailDesign.add_command(label='垂尾根弦长和梢弦长', command=self.callback_air_chuiwei_xianchang)
        airDesign.add_cascade(label="垂尾", menu=verticalTailDesign)

        # 重量分析子菜单
        WeightDesign = Menu(emenu, tearoff=False)
        WeightDesign.add_command(label='机翼重量估算', command=self.callback_weight_wing)
        WeightDesign.add_command(label='机身重量估算', command=self.callback_weight_fuselage)
        WeightDesign.add_command(label='尾翼重量估算', command=self.callback_weight_tail)
        WeightDesign.add_command(label='起落架重量估算', command=self.callback_weight_landing)

        # 结构设计子菜单
        structDesign = Menu(emenu, tearoff=False)

        overloadDesign = Menu(structDesign, tearoff=False)
        overloadDesign.add_command(label='俯冲拉起', command=self.callback_FCLQ)
        overloadDesign.add_command(label='等速盘旋', command=self.callback_DSPX)
        overloadDesign.add_command(label='遭遇垂直突风', command=self.callback_CZTF)
        structDesign.add_cascade(label='过载系数', menu=overloadDesign)

        strucCheck = Menu(structDesign, tearoff=False)
        strucCheck.add_command(label='翼梁强度校核', command=self.callback_QD)
        strucCheck.add_command(label='翼梁稳定性校核', command=self.callback_WDX)
        structDesign.add_cascade(label='结构校核', menu=strucCheck)

        # 性能设计子菜单
        performanceDesign = Menu(emenu, tearoff=False)

        #升阻性能
        airAnalysis = Menu(performanceDesign, tearoff=False)
        airAnalysis.add_command(label='机翼升力线斜率', command=self.callback_SLXXL)
        airAnalysis.add_command(label='奥斯瓦尔德效率因子', command=self.callback_Oswald)
        airAnalysis.add_command(label='巡航升力系数', command=self.callback_PF)
        airAnalysis.add_command(label='诱导阻力系数', command=self.callback_DI)
        airAnalysis.add_command(label='最大升阻比', command=self.callback_E_max)
        performanceDesign.add_cascade(label='升阻性能', menu=airAnalysis)

        #起降性能
        TAKEOFF_LAND = Menu(performanceDesign, tearoff=False)
        TAKEOFF_LAND.add_command(label='起飞离地速度', command=self.callback_TAKEOFF_velocity)
        TAKEOFF_LAND.add_command(label='起飞滑跑距离', command=self.callback_TAKEOFF_distance)
        TAKEOFF_LAND.add_command(label='进场速度', command=self.callback_LAND_velocity)
        performanceDesign.add_cascade(label='起降性能', menu=TAKEOFF_LAND)

        #续航性能
        Continue = Menu(performanceDesign, tearoff=False)
        Continue.add_command(label='等速巡航航程(喷气式)', command=self.callback_Cruise_distance)
        Continue.add_command(label='航时(喷气式)', command=self.callback_Cruise_time)
        performanceDesign.add_cascade(label='续航性能', menu=Continue)

        #盘旋性能
        HOVERING = Menu(performanceDesign, tearoff=False)
        HOVERING.add_command(label='性能指标(水平协调盘旋)', command=self.callback_hover)
        performanceDesign.add_cascade(label='盘旋性能', menu=HOVERING)

        # 稳定性设计子菜单
        stabilityDesign = Menu(emenu, tearoff=False)
        stabilityDesign.add_command(label='纵向静稳定裕度', command=self.callback_H_n)

        # 大气参数子菜单
        costDesign = Menu(emenu, tearoff=False)
        costDesign.add_command(label='大气参数', command=self.callback_air)
        costDesign.add_command(label='翼型升力', command=self.callback_lift_airfoil)
        costDesign.add_command(label='马赫角', command=self.callback_mahe_angel)
        costDesign.add_command(label='平板层流边界层厚度', command=self.callback_boundary_layer)
        costDesign.add_command(label='伯努利方程(皮托管测空速)', command=self.callback_bonuli)
        costDesign.add_command(label='线化超音速理论(小迎角平板翼型)', command=self.callback_linear)

        # for item in ['气动设计', '结构设计', '性能分析', '稳定性分析', '费用分析']:
        #     emenu.add_command(label=item, command = self.callback)

        # 在 menubar 上设置子菜单名，并关联一系列子菜单
        menubar.add_cascade(label="文件", menu=fmenu)
        menubar.add_cascade(label="总体设计", menu=mainDesign)
        menubar.add_cascade(label="气动设计", menu=airDesign)
        menubar.add_cascade(label="重量估算", menu=WeightDesign)
        menubar.add_cascade(label="结构设计", menu=structDesign)
        menubar.add_cascade(label="性能分析", menu=performanceDesign)
        menubar.add_cascade(label="稳定性分析", menu=stabilityDesign)
        menubar.add_cascade(label="空气动力学", menu=costDesign)

        # 显示菜单
        self.root.config(menu=menubar)
        self.inputFrame = Frame(self.root, height=200, width=1000, highlightthickness=2, relief=SUNKEN,
                                highlightbackground='#000000')
        # 输入标签
        self.input_tip_label = Label(self.inputFrame, text="输入", font='50px')
        self.input_tip_label.grid(row=0, column=0, sticky='n')
        self.inputLabel = [None, None, None, None, None, None, None, None, None, None]
        self.inputLabel[0] = Label(self.inputFrame, text="空机重量", font='50px')
        self.inputLabel[1] = Label(self.inputFrame, text="乘员重量", font='50px')
        self.inputLabel[2] = Label(self.inputFrame, text="装载重量", font='50px')
        self.inputLabel[3] = Label(self.inputFrame, text="燃油重量", font='50px')
        self.show_input_lable()

        # 输入文本框
        self.inputText = [None, None, None, None, None, None, None, None, None, None]
        self.inputText[0] = Entry(self.inputFrame, width=75, font='50px')
        self.inputText[1] = Entry(self.inputFrame, width=75, font='50px')
        self.inputText[2] = Entry(self.inputFrame, width=75, font='50px')
        self.inputText[3] = Entry(self.inputFrame, width=75, font='50px')
        self.show_input_text()

        # 输入单位
        self.inputUnit = [None, None, None, None, None, None, None,None, None, None]
        self.inputUnit[0] = Label(self.inputFrame, text='kg', font='50px')
        self.inputUnit[1] = Label(self.inputFrame, text='kg', font='50px')
        self.inputUnit[2] = Label(self.inputFrame, text='kg', font='50px')
        self.inputUnit[3] = Label(self.inputFrame, text='kg', font='50px')
        self.show_input_unit()

        # 总输入内容绘制
        self.inputFrame.grid(row=0, column=0, padx=40, pady=10, ipady=10)

        # 计算按钮
        self.get_button = Button(self.root, text="计算", bg="#EEEEEE", width=10,
                                 command=self.get_result)  # 调用内部方法  加()为直接调用
        self.get_button.grid(row=1, column=0)

        # 输出标签
        self.outputFrame = Frame(self.root, height=200, width=1000, highlightthickness=2, relief=SUNKEN,
                                 highlightbackground='#000000')
        self.output_tip_label = Label(self.outputFrame, text="输出", font='50px')
        self.output_tip_label.grid(row=0, column=0, sticky='n')
        self.outputLabel = [None, None, None, None, None, None, None, None,None, None, None]
        self.outputLabel[0] = Label(self.outputFrame, text="起飞重量", font='50px')
        self.show_output_lable()

        # 输出文本
        self.outputText = [None, None, None, None, None, None, None, None,None, None, None]
        self.outputText[0] = Text(self.outputFrame, width=75, height=1, font='50px', bg='#FFFFFF', highlightthickness=1,
                                  relief='flat')
        # self.outputText[0].insert("insert", "50")
        self.show_output_text()

        # 输出单位
        self.outputUnit = [None, None, None, None, None, None, None, None,None, None, None]
        self.outputUnit[0] = Label(self.outputFrame, text='kg', font='50px')
        self.show_output_unit()

        # 输出总绘制
        self.outputFrame.grid(row=2, column=0, padx=40, pady=10, ipady=10)

        # 日志
        self.logFrame = Frame(self.root, height=200, width=1000, highlightthickness=2, relief=SUNKEN,
                              highlightbackground='#000000')
        self.log_tip_label = Label(self.logFrame, text="日志", font='50px')
        self.log_tip_label.grid(row=0, column=0, sticky='n')

        # 输出文本
        self.logText = Text(self.logFrame, height=10, font='50px', bg='#EEEEEE', highlightthickness=0,
                            relief='flat')
        self.logText.grid()

        # 输出总绘制
        self.logFrame.grid(row=3, column=0, padx=40, pady=10, ipady=10)

        # 计算按钮
        self.lh_button = Button(self.root, text="上一条", bg="#EEEEEE", width=10,
                                command=self.last_history)  # 调用内部方法  加()为直接调用
        self.lh_button.grid(row=4, column=0, pady=5)
        self.nh_button = Button(self.root, text="下一条", bg="#EEEEEE", width=10,
                                command=self.next_history)  # 调用内部方法  加()为直接调用
        self.nh_button.grid(row=5, column=0, pady=5)

    # 功能函数
    def str_trans_to_md5(self):
        src = self.init_data_Text.get(1.0, END).strip().replace("\n", "").encode()
        # print("src =",src)
        if src:
            try:
                myMd5 = hashlib.md5()
                myMd5.update(src)
                myMd5_Digest = myMd5.hexdigest()
                # print(myMd5_Digest)
                # 输出到界面
                self.result_data_Text.delete(1.0, END)
                self.result_data_Text.insert(1.0, myMd5_Digest)
                self.write_log_to_Text("INFO:str_trans_to_md5 success")
            except:
                self.result_data_Text.delete(1.0, END)
                self.result_data_Text.insert(1.0, "字符串转MD5失败")
        else:
            self.write_log_to_Text("ERROR:str_trans_to_md5 failed")

    # 获取当前时间
    def get_current_time(self):
        current_time = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time()))
        return current_time

    # 日志动态打印
    def write_log_to_Text(self, logmsg):
        global LOG_LINE_NUM
        current_time = self.get_current_time()
        logmsg_in = str(current_time) + " " + str(logmsg) + "\n"  # 换行
        if LOG_LINE_NUM <= 7:
            self.log_data_Text.insert(END, logmsg_in)
            LOG_LINE_NUM = LOG_LINE_NUM + 1
        else:
            self.log_data_Text.delete(1.0, 2.0)
            self.log_data_Text.insert(END, logmsg_in)

def gui_start():
    init_window = Tk()  # 实例化出一个父窗口
    ZMJ_PORTAL = MY_GUI(init_window)
    # 设置根窗口默认属性
    ZMJ_PORTAL.set_init_window()

    init_window.mainloop()  # 父窗口进入事件循环，可以理解为保持窗口运行，否则界面不展示


gui_start()
