from cProfile import label
import scipy.stats as stats
from sympy import *
import numpy as np
from rich import print
import math
import xlwt

'''
    求解二项式分布可靠度单侧置信下限
'''

def beize_pulate():
    '''贝泽普拉特(其它场合)'''

    n = 4
    F = 1
    r = 0.50

    ur = stats.norm.ppf(r)
    R = symbols('R', real = True)
    d = (F + 0.5 + 1/6) - (n+1/3)*(1-R) + 0.02*(R/(F+1) - (1-R)/(n-F) + (R-0.5)/(n+1))
    h = ( (d/-(F+0.5-n*(1-R))) * sqrt((2/(1+(1/(6*n)))) * ((F+0.5)*log((F+0.5)/(n*(1-R))) + (n-F-0.5)*log((n-F-0.5)/(n*R)))) ) + ur

    EQ = h 
    a = nsolve(EQ, 0.5, solver='newton', verify=False)

    import numpy as np
    import matplotlib.pyplot as plt
    xx = np.linspace(-10.0, 10.0, 1000)
    func = lambdify(R, h, 'numpy')
    yy = [ func(x) for x in xx]
    plt.plot(xx,yy)

    plt.grid(True)
    plt.show()

def draw(R, h):

    import numpy as np
    import matplotlib.pyplot as plt
    xx = np.linspace(-10.0, 10.0, 1000)
    func = lambdify(R, h, 'numpy')
    yy = [ func(x) for x in xx]
    plt.plot(xx,yy)
    plt.grid(True)
    plt.show()

def log_gamma(r, n, F):
    '''对数伽马近似(F = 1,2,3)'''
    
    c = math.log((n+1)/(n-F)) / math.log((n+2)/(n-F+1))
    μ = (3-c) / (2*(c-1) - 0.355*(c-1)**3)
    Z = math.log((n+1)/(n-F)) / (math.log((μ+1)/μ))
    chi2 = stats.chi2.ppf(r, df=2*Z)

    return math.exp(-chi2/(2*μ))

def beize_pulate2(r, n, F):
    '''贝泽普拉特'''

    ur = stats.norm.ppf(r)
    R = symbols('R', real = True)
    d = (F + 0.5 + 1/6) - (n+1/3)*(1-R) + 0.02*(R/(F+1) - (1-R)/(n-F) + (R-0.5)/(n+1))
    h = ( (d/-(F+0.5-n*(1-R))) * sqrt((2/(1+(1/(6*n)))) * ((F+0.5)*log((F+0.5)/(n*(1-R))) + (n-F-0.5)*log((n-F-0.5)/(n*R)))) ) + ur
    try:
        a = nsolve(h, 0.4, solver='newton', verify=False)
        return a
    except:
        return 0

def find_binomial_distribution_one_side_confidence_lower_limit_table(r, n, F):
    # 约束条件：n > 0, F >= 0，r ={0.50, 0.60, 0.70, 0.80, 0.90, 0.95, 0.99} 
    if F == 0:
        return pow(1-r, 1/n)
    elif F >= n:
        value = 0
    elif ((2*F)+1)  == n:
        value = 0.5
    elif F == n-1:
        return 1-pow(r, 1/n)
    elif F in {1,2,3}:
        value = log_gamma(r, n, F)
    else:
        value = beize_pulate2(r, n, F)

    if type(value) == Add:
        value = re(value)

    print("n: {}, F: {} ----> RL: {:.7f}".format(n, F, value))

    return value


# 测试F=0情况
def test_F_is_zero_case():

    r, F = 0.50, 0.0
    print(f"============ 置信水平: {r} ============ ")
    for n in range(1, 1001):
        find_binomial_distribution_one_side_confidence_lower_limit_table(r, n, F)

def test_F_is_one_case():

    r, F = 0.50, 1
    print(f"============ 置信水平: {r} ============ ")
    for n in range(1, 1001):
        find_binomial_distribution_one_side_confidence_lower_limit_table(r, n, F)

def test_F_is_tow_case():

    r, F = 0.50, 2
    print(f"============ 置信水平: {r} ============ ")
    for n in range(1, 1001):
        find_binomial_distribution_one_side_confidence_lower_limit_table(r, n, F)


def test_F_is_three_case():

    r, F = 0.50, 3
    print(f"============ 置信水平: {r} ============ ")
    for n in range(1, 1001):
        find_binomial_distribution_one_side_confidence_lower_limit_table(r, n, F)
        
def test_other_case():

    r, F = 0.50, 20
    print(f"============ 置信水平: {r} ============ ")
    for n in np.arange(1, 70, 1):
        find_binomial_distribution_one_side_confidence_lower_limit_table(r, n, F)

    for n in np.arange(70, 100, 5):
        find_binomial_distribution_one_side_confidence_lower_limit_table(r, n, F)

    for n in np.arange(100, 140, 10):
        find_binomial_distribution_one_side_confidence_lower_limit_table(r, n, F)

    for n in np.arange(140, 200, 20):
        find_binomial_distribution_one_side_confidence_lower_limit_table(r, n, F)

    for n in np.arange(200, 1001, 50):
        find_binomial_distribution_one_side_confidence_lower_limit_table(r, n, F)

# test_other_case()

def do(r, n, F):
    return find_binomial_distribution_one_side_confidence_lower_limit_table(r, n, F)

def gen_data_txt_file():

    with open('./data.txt', 'wb') as f:
        for r in {0.50, 0.60, 0.70, 0.80, 0.90, 0.95, 0.99}:
            # str2 = "============ 置信水平: " + str(r) + "\n"
            # f.write(str2.encode())

            for F in range(0, 21):
                for n in np.arange(1, 70, 1):
                    str3 = "r: {:.2f}, n: {}, F: {}, RL: {:.7f} \n".format(r, n, F, do(r, n, F))
                    f.write(str3.encode())

                # for n in np.arange(70, 100, 5):
                #     str3 = "r: {:.2f}, n: {}, F: {}, RL: {:.7f} \n".format(r, n, F, do(r, n, F))
                #     f.write(str3.encode())

                # for n in np.arange(100, 140, 10):
                #     str3 = "r: {:.2f}, n: {}, F: {}, RL: {:.7f} \n".format(r, n, F, do(r, n, F))
                #     f.write(str3.encode())

                # for n in np.arange(140, 200, 20):
                #     str3 = "r: {:.2f}, n: {}, F: {}, RL: {:.7f} \n".format(r, n, F, do(r, n, F))
                #     f.write(str3.encode())

                # for n in np.arange(200, 1001, 50):
                #     str3 = "r: {:.2f}, n: {}, F: {}, RL: {:.7f} \n".format(r, n, F, do(r, n, F))
                #     f.write(str3.encode())

        # wb.save('binomial_distribution_one_side_confidence_lower_limit_table.xls').

# gen_data_txt_file()

def gen_data_excel_file():

    wb = xlwt.Workbook(encoding = 'utf-8')
    worksheet = wb.add_sheet('xxxx', cell_overwrite_ok=True)
    alignment = xlwt.Alignment()
    alignment.horz = xlwt.Alignment.HORZ_CENTER
    style =xlwt.XFStyle()
    style.alignment = alignment
    
    style2 =xlwt.XFStyle()
    font = xlwt.Font() # 为样式创建字体
    font.bold = True
    # font.size = 28
    style2.font = font

    worksheet.col(0).width = 256 * 5

    row_idx = 0
    for r in [0.50, 0.60, 0.70, 0.80, 0.90, 0.95, 0.99]:
    # for r in [0.50]:

        worksheet.write(row_idx, 0, label = 'γ = {:.2f}'.format(r), style = style2)
        row_idx += 1

        for i in range(0, 21):
            worksheet.col(i+1).width = 256 * 12
            worksheet.write(row_idx, i+1, label = i, style = style)
        row_idx += 1

        # 1 ~ 70
        for n in np.arange(1, 70, 1):
            col_idx = 0
            for F in range(0, 21):
                if col_idx == 0: 
                    worksheet.write(row_idx, col_idx, label = str(n), style = style)

                worksheet.write(row_idx, col_idx+1, label = '{:.7f}'.format(do(r, n, F)), style = style)
                col_idx += 1
            row_idx += 1

        # 70 ~ 100
        for n in np.arange(70, 100, 5):
            col_idx = 0
            for F in range(0, 21):
                if col_idx == 0: 
                    worksheet.write(row_idx, col_idx, label = str(n), style = style)

                worksheet.write(row_idx, col_idx+1, label = '{:.7f}'.format(do(r, n, F)), style = style)
                col_idx += 1
            row_idx += 1

        # 100 ~ 140
        for n in np.arange(100, 140, 10):
            col_idx = 0
            for F in range(0, 21):
                if col_idx == 0: 
                    worksheet.write(row_idx, col_idx, label = str(n), style = style)

                worksheet.write(row_idx, col_idx+1, label = '{:.7f}'.format(do(r, n, F)), style = style)
                col_idx += 1
            row_idx += 1

        # 140 ~ 200
        for n in np.arange(140, 200, 20):
            col_idx = 0
            for F in range(0, 21):
                if col_idx == 0: 
                    worksheet.write(row_idx, col_idx, label = str(n), style = style)

                worksheet.write(row_idx, col_idx+1, label = '{:.7f}'.format(do(r, n, F)), style = style)
                col_idx += 1
            row_idx += 1
        # 200 ~ 1001
        for n in np.arange(200, 1001, 50):
            col_idx = 0
            for F in range(0, 21):
                if col_idx == 0: 
                    worksheet.write(row_idx, col_idx, label = str(n), style = style)

                worksheet.write(row_idx, col_idx+1, label = '{:.7f}'.format(do(r, n, F)), style = style)
                col_idx += 1
            row_idx += 1

    wb.save('二项分布可靠度单侧置信下限数表.xls')

gen_data_excel_file()
