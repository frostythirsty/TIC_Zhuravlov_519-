import matplotlib.pyplot as plt
import numpy
import numpy as np
import scipy
#9
n = 500
Fs = 1000
F_max = 19
F_filter = 26
random = numpy.random.normal(0, 10, n)
x = numpy.arange(n)/Fs
w = F_max/(Fs/2)
parameter = scipy.signal.butter(3, w, 'low', output='sos')
y = scipy.signal.sosfiltfilt(parameter, random)
fig,ax = plt.subplots(figsize=(21/2.54, 14/2.54))
ax.plot(x, y, linewidth=1)
ax.set_xlabel("Час(секунди)", fontsize=14)
ax.set_ylabel("Амплітуда сигналу", fontsize=14)
plt.title("Сигнал з максимальною частотою F_max = 31 Гц", fontsize=14)
fig.savefig('./figures/'+'Сигнал з максимальною частотою F_max = 31 Гц'+'.png', dpi=600)
plt.close(fig)
y_s = numpy.abs(scipy.fft.fftshift(scipy.fft.fft(y)))
x_s = scipy.fft.fftshift(scipy.fft.fftfreq(n, 1/n))
fig,ax = plt.subplots(figsize=(21/2.54, 14/2.54))
ax.plot(x_s, y_s, linewidth=1)
ax.set_xlabel("Частота (Герци)", fontsize=14)
ax.set_ylabel("Амплітуда спектру", fontsize=14)
plt.title("Спектр сигналу", fontsize=14)
fig.savefig('./figures/'+'Спектр сигналу'+'.png', dpi=600)
plt.close(fig)
discrete_signals = []
discrete_spectrums = []
signal_after_filers = []
dispersion = []
variance_dif = []
for Dt in [2, 4, 8, 16]:
    discrete_signal = numpy.zeros(n)
    for i in range(0, round(n/Dt)):
        discrete_signal[i * Dt] = y[i * Dt]
    discrete_signals += [list(discrete_signal)]
    y_s = numpy.abs(scipy.fft.fftshift(scipy.fft.fft(discrete_signal)))
    discrete_spectrums += [list(y_s)]
    w = F_filter/(Fs/2)
    parameters_filter = scipy.signal.butter(3, w, 'low', output='sos')
    discrete_signal_after_filers = scipy.signal.sosfiltfilt(parameters_filter, discrete_signal)
    signal_after_filers += [list(discrete_signal_after_filers)]
    E1 = discrete_signal_after_filers - y
    dispersion += [numpy.var(E1)]
    variance_dif += [numpy.var(y) / numpy.var(E1)]
fig, ax = plt.subplots(2, 2, figsize=(21/2.54, 14/2.54))
s = 0
for i in range(0, 2):
    for j in range(0, 2):
        ax[i][j].plot(x, discrete_signals[s], linewidth=1)
        s += 1
fig.supxlabel('Час(секунди)', fontsize=14)
fig.supylabel('Амплітуда сигналу', fontsize=14)
fig.suptitle('Сигнал з кроком дискретизації Dt = (2, 4, 8, 16)', fontsize=14)
fig.savefig('./figures/'+'Сигнал з кроком дискретизації'+'.png', dpi=600)
plt.close(fig)
fig, ax = plt.subplots(2, 2, figsize=(21/2.54, 14/2.54))
s = 0
for i in range(0, 2):
    for j in range(0, 2):
        ax[i][j].plot(x_s, discrete_spectrums[s], linewidth=1)
        s += 1
fig.supxlabel('Частота(Гц)', fontsize=14)
fig.supylabel('Амплітуда спектру', fontsize=14)
fig.suptitle('Спектри сигналів з кроком дискретизації Dt = (2, 4, 8, 16)', fontsize=14)
fig.savefig('./figures/' + 'Спектр сигналів з кроком дискретизації' + '.png', dpi=600)
plt.close(fig)
fig, ax = plt.subplots(2, 2, figsize=(21/2.54, 14/2.54))
s = 0
for i in range(0, 2):
    for j in range(0, 2):
        ax[i][j].plot(x, signal_after_filers[s], linewidth=1)
        s += 1
fig.supxlabel('Час(секунди)', fontsize=14)
fig.supylabel('Амплітуда сигналу', fontsize=14)
fig.suptitle('Відновлені аналогові сигнали з кроком дискретизації Dt = (2, 4, 8, 16)', fontsize=14)
fig.savefig('./figures/' + 'Відновлені аналогові сигнали з кроком дискретизації' + '.png', dpi=600)
plt.close(fig)
fig,ax = plt.subplots(figsize=(21/2.54, 14/2.54))
ax.plot([2, 4, 8, 16], dispersion, linewidth=1)
fig.supxlabel('Крок дискретизації', fontsize=14)
fig.supylabel('Дисперсія', fontsize=14)
fig.suptitle('Залежність дисперсії від кроку дискретизації', fontsize=14)
fig.savefig('./figures/' + 'Залежність дисперсії від кроку дискретизації' + '.png', dpi=600)
plt.close(fig)
fig,ax = plt.subplots(figsize=(21/2.54, 14/2.54))
ax.plot([2, 4, 8, 16], variance_dif, linewidth=1)
fig.supxlabel('Крок дискретизації', fontsize=14)
fig.supylabel('ССШ', fontsize=14)
fig.suptitle('Залежність співвідношення сигнал-шум від кроку дискретизації', fontsize=14)
fig.savefig('./figures/' + 'Залежність співвідношення сигнал-шум від кроку дискретизації' + '.png', dpi=600)
plt.close(fig)
quantize_signal = []
dispers = []
signal_to_noise_ratio = []
quantize_signals = []
for M in [4, 16, 64, 256]:
    bits = []
    bit_signal = []
    delta = (numpy.max(y) - numpy.min(y)) / (M - 1) # крок квантування
    quantize_signal = delta * np.round(y / delta) # цифровій сигнал
    quantize_levels = numpy.arange(numpy.min(quantize_signal), numpy.max(quantize_signal)+1, delta) # список дискретних відліків амплітуди сигналу з кроком delta
    quantize_signals += [quantize_signal]
    quantize_bit = numpy.arange(0, M)
    quantize_bit = [format(bits, '0' + str(int(numpy.log(M) / numpy.log(2))) + 'b') for bits in quantize_bit] # перетворення десяткових значень у біти
    quantize_table = numpy.c_[quantize_levels[:M], quantize_bit[:M]]
    fig, ax = plt.subplots(figsize=(14 / 2.54, M / 2.54))
    table = ax.table(cellText=quantize_table, colLabels=['Значення сигналу','Кодова послідовність'], loc='center')
    table.set_fontsize(14)
    table.scale(1, 2)
    ax.axis('off')
    fig.savefig('./figures/' + 'Таблиця квантування для ' + str(M) + ' рівнів' + '.png', dpi=600)
    plt.close(fig)
    for signal_value in quantize_signal:
        for index, value in enumerate(quantize_levels[:M]):
            if numpy.round(numpy.abs(signal_value - value), 0) == 0:
                bits.append(quantize_bit[index])
                break

    bits = [int(item) for item in list(''.join(bits))]
    fig, ax = plt.subplots(figsize=(21 / 2.54, 14 / 2.54))
    x_b = numpy.arange(0, len(bits))
    ax.step(x_b, bits, linewidth=0.1)
    fig.supxlabel('Біти', fontsize=14)
    fig.supylabel('Амплітуда сигналу', fontsize=14)
    fig.suptitle('Кодова послідовність сигналу при кількості рівнів квантування 256', fontsize=14)
    fig.savefig('./figures/' + 'Кодова послідовність сигналу при кількості рівнів квантування ' + str(M) + '.png', dpi=600)
    plt.close(fig)
    E2 = quantize_signal - y
    dispers += [numpy.var(E2)]
    signal_to_noise_ratio += [numpy.var(y) / numpy.var(E2)]
fig, ax = plt.subplots(2, 2, figsize=(21/2.54, 14/2.54))
s = 0
for i in range(0, 2):
    for j in range(0, 2):
        ax[i][j].plot(x, quantize_signals[s], linewidth=1)
        s += 1
fig.supxlabel('Час(секунди)', fontsize=14)
fig.supylabel('Амплітуда сигналу', fontsize=14)
fig.suptitle('Цифрові сигнали з рівнями квантування (4, 16, 64, 256)', fontsize=14)
fig.savefig('./figures/' + 'Цифрові сигнали з рівнями квантування (4, 16, 64, 256)' + '.png', dpi=600)
plt.close(fig)
fig,ax = plt.subplots(figsize=(21/2.54, 14/2.54))
ax.plot([4, 16, 64, 256], dispers, linewidth=1)
fig.supxlabel('Кількість рівнів квантування', fontsize=14)
fig.supylabel('Дисперсія', fontsize=14)
fig.suptitle('Залежність дисперсії від кількості рівніів квантування', fontsize=14)
fig.savefig('./figures/' + 'Залежність дисперсії від кількості рівніів квантування' + '.png', dpi=600)
plt.close(fig)
fig,ax = plt.subplots(figsize=(21/2.54, 14/2.54))
ax.plot([4, 16, 64, 256], signal_to_noise_ratio, linewidth=1)
fig.supxlabel('Кількість рівнів квантування', fontsize=14)
fig.supylabel('ССШ', fontsize=14)
fig.suptitle('Залежність співвідношення сигнал-шум від кількості рівніів квантування', fontsize=14)
fig.savefig('./figures/' + 'Залежність співвідношення сигнал-шум від кількості рівніів квантування' + '.png', dpi=600)
plt.close(fig)