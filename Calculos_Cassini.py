import numpy as np
import matplotlib.pyplot as plt
import os

mu_Sol = 0.0002959122

mu_planetas = {
    'Venus': 0.0000024478383,
    'Terra': 0.0000030034896,
    'Jupiter': 0.0009547919,
    'Saturno': 0.0002858850
}

eventos_jd = {
    'Lança_Terra': 2450737.5,
    'Venus1':      2450929.5,
    'Venus2':      2451353.5,
    'Terra':       2451408.5,
    'Jupiter':     2451908.5,
    'Saturno':     2453187.5
}

cores_fases = ["blue", 'orange', 'green', 'red', 'purple']

def rotacionar_90_horario(x, y):
    x_rot = y
    y_rot = -x
    return x_rot, y_rot

class Orbita:
    def __init__(self):
        self.a = 0
        self.e = 0
        self.p = 0
        self.omega = 0
        self.rvec = np.zeros(3)
        self.vvec = np.zeros(3)
        self.e_vec = np.zeros(3)

    def vet_estado(self, mu_Sol, rvec, vvec):
        self.rvec = rvec
        self.vvec = vvec
        r = np.linalg.norm(rvec)
        v = np.linalg.norm(vvec)
        self.a = 1 / ((2 / r) - (v**2 / mu_Sol))
        hvec = np.cross(rvec, vvec)
        termo1 = (v**2 - (mu_Sol / r)) * rvec
        termo2 = np.dot(rvec, vvec) * vvec
        self.e_vec = (1 / mu_Sol) * (termo1 - termo2)
        self.e = np.linalg.norm(self.e_vec)
        self.p = np.linalg.norm(hvec)**2 / mu_Sol
        self.omega = np.arctan2(self.e_vec[1], self.e_vec[0])

    def raio_orbita(self, f):
        return (self.p) / (1 + self.e * np.cos(f))

    def pos_orb(self, f):
        r = self.raio_orbita(f)
        xpq = r * np.cos(self.omega + f)
        ypq = r * np.sin(self.omega + f)
        return xpq, ypq

    def obter_anomalia_verdadeira(self, r_vec_alvo):
        theta_alvo = np.arctan2(r_vec_alvo[1], r_vec_alvo[0])
        f = theta_alvo - self.omega
        f = (f + np.pi) % (2 * np.pi) - np.pi
        return f

    def flyby(self, r_planeta, v_planeta, mu_planeta):
        delta_r = self.rvec - r_planeta
        rperiastro = np.linalg.norm(delta_r)
        v_infty_vec = self.vvec - v_planeta
        v_infty = np.linalg.norm(v_infty_vec)
        e_hip = 1 + (rperiastro * v_infty**2) / mu_planeta
        delta = 2*np.arcsin(1/e_hip)
        return rperiastro, v_infty, e_hip, np.degrees(delta)

    def ganho_flyby(self, v_infty, delta_deg):
        delta_rad = np.radians(delta_deg)
        return 2 * v_infty * np.sin(delta_rad / 2)

def ler_horizons_completo(caminho_arquivo):
    dados = []
    lendo = False
    try:
        with open(caminho_arquivo, 'r') as f:
            for line in f:
                line = line.strip()
                if line == "$$SOE":
                    lendo = True
                    continue
                if line == "$$EOE":
                    lendo = False
                    continue
                if lendo:
                    try:
                        parts = line.split(',')
                        jd = float(parts[0])
                        r = np.array([float(parts[2]), float(parts[3]), float(parts[4])])
                        v = np.array([float(parts[5]), float(parts[6]), float(parts[7])])
                        dados.append({'jd': jd, 'r': r, 'v': v})
                    except ValueError:
                        continue
    except FileNotFoundError:
        return []
    return dados

# Carregamento de dados
arquivos_planetas = {'Venus': 'VetoresVenus.txt', 'Terra': 'VetoresTerra.txt', 'Jupiter': 'VetoresJupiter.txt', 'Saturno': 'VetoresSaturno.txt'}
dados_planetas = {}
for nome, arquivo in arquivos_planetas.items():
    dados = ler_horizons_completo(os.path.join("Dados", arquivo))
    if dados:
        dados_planetas[nome] = dados

dados_cassini = ler_horizons_completo(os.path.join("Dados", "VetoresCassini.txt"))

# Plotagem
plt.figure(figsize=(12, 12))
plt.title("Trajetória da Cassini")
plt.xlabel("Distância [AU]")
plt.ylabel("Distância [AU]")
plt.axis('equal')
plt.grid(True, linestyle='--', alpha=0.3)
plt.plot(0, 0, 'yo', markersize=12, label='Sol')

fases = [
    ('Terra -> Venus 1', eventos_jd['Lança_Terra'], eventos_jd['Venus1']),
    ('Venus 1 -> Venus 2', eventos_jd['Venus1'], eventos_jd['Venus2']),
    ('Venus 2 -> Terra', eventos_jd['Venus2'], eventos_jd['Terra']),
    ('Terra -> Jupiter', eventos_jd['Terra'], eventos_jd['Jupiter']),
    ('Jupiter -> Saturno', eventos_jd['Jupiter'], eventos_jd['Saturno'])
]

distancia_total_AU = 0

if dados_cassini:
    jds = np.array([d['jd'] for d in dados_cassini])

    for i, (nome_fase, jd_inicio, jd_fim) in enumerate(fases):
        cor = cores_fases[i % len(cores_fases)]

        indices = np.where((jds >= jd_inicio) & (jds <= jd_fim))[0]

        if len(indices) > 0:
            # Dados para plotar
            segmento_r = np.array([dados_cassini[k]['r'] for k in indices])
            x_real, y_real = segmento_r[:, 0], segmento_r[:, 1]
            x_plot, y_plot = rotacionar_90_horario(x_real, y_real)

            plt.plot(x_plot, y_plot, '-', color=cor, linewidth=1.5, label=nome_fase)
            plt.plot(x_plot[0], y_plot[0], 'o', color=cor, markersize=4)

            ponto_inicio = dados_cassini[indices[0]]
            ponto_fim = dados_cassini[indices[-1]]

            orb = Orbita()
            orb.vet_estado(mu_Sol, ponto_inicio['r'], ponto_inicio['v'])

            f_start = orb.obter_anomalia_verdadeira(ponto_inicio['r'])
            f_end = orb.obter_anomalia_verdadeira(ponto_fim['r'])
            if f_end < f_start:
                f_end += 2 * np.pi

            # Plotagem da órbita completa tracejada
            f_full = np.linspace(0, 2 * np.pi, 300)
            x_full, y_full = orb.pos_orb(f_full)
            x_full_rot, y_full_rot = rotacionar_90_horario(x_full, y_full)
            plt.plot(x_full_rot, y_full_rot, ':', color=cor, alpha=0.5, linewidth=1)

# Adiciona a chegada em Saturno
if 'Saturno' in dados_planetas:
    jd_chegada = eventos_jd['Saturno']
    jds_s = np.array([d['jd'] for d in dados_planetas['Saturno']])
    idx_s = (np.abs(jds_s - jd_chegada)).argmin()
    pt_saturno = dados_planetas['Saturno'][idx_s]

    orb_sat = Orbita()
    orb_sat.vet_estado(mu_Sol, pt_saturno['r'], pt_saturno['v'])


    f_sat = np.linspace(0, 2 * np.pi, 300)
    x_sat, y_sat = orb_sat.pos_orb(f_sat)
    x_s_rot, y_s_rot = rotacionar_90_horario(x_sat, y_sat)

    plt.plot(x_s_rot, y_s_rot, '--', color='brown', alpha=0.6, linewidth=1, label='Órbita Saturno')

    xs_real, ys_real = rotacionar_90_horario(pt_saturno['r'][0], pt_saturno['r'][1])
    plt.plot(xs_real, ys_real, 'o', color='brown', markersize=5, label='Saturno')

    theta_anel = np.linspace(0, 2 * np.pi, 100)

    # Raio visual do anel
    raio_anel = 0.25

    x_anel = xs_real + raio_anel * np.cos(theta_anel)
    y_anel = ys_real + (raio_anel) * np.sin(theta_anel) # Multiplico por 0.3 para "achatar" (efeito 3D)

    # Plotar o anel
    plt.plot(x_anel, y_anel, '-', color='peru', alpha=0.8, linewidth=1.5)

plt.legend(loc='upper left', fontsize='small', framealpha=0.9, bbox_to_anchor=(1, 1))
plt.tight_layout()
plt.xlim(-2,10)
plt.ylim(-4,6)
plt.savefig("./Outputs/trajetoria_completa_cassini.png", dpi=300)
plt.show()

# Calculo Flybys
ocorre_flyby = [('Venus1', 'Venus'), ('Venus2', 'Venus'), ('Terra', 'Terra'), ('Jupiter', 'Jupiter')]

with open("./Outputs/resultados_flyby.txt", "w") as f:
    f.write("Resultados dos Flybys:\n\n")

for flyby, planeta in ocorre_flyby:
    jd_flyby = eventos_jd[flyby]

    idx_flyby = (np.abs(jds - jd_flyby)).argmin()
    ponto_cassini = dados_cassini[idx_flyby]

    jds_planeta = np.array([d['jd'] for d in dados_planetas[planeta]])
    idx_p = (np.abs(jds_planeta - jd_flyby)).argmin()
    ponto_planeta = dados_planetas[planeta][idx_p]

    orb_calc = Orbita()
    orb_calc.vet_estado(mu_Sol, ponto_cassini['r'], ponto_cassini['v'])

    r_p, v_inf, e_h, delta = orb_calc.flyby(ponto_planeta['r'], ponto_planeta['v'], mu_planetas[planeta])
    ganho = orb_calc.ganho_flyby(v_inf, delta)
    ganho_km = ganho * 1731.45683

    with open("./Outputs/resultados_flyby.txt", "a") as f:
        f.write(f"Flyby em {planeta} ({flyby}):\n")
        f.write(f"  - V_infty: {v_inf:.6f} AU/dia\n")
        f.write(f"  - Excentricidade: {e_h:.6f}\n")
        f.write(f"  - Deflexao: {delta:.6f} graus\n")
        f.write(f"  - Ganho V: {ganho_km:.6f} km/s\n\n")
