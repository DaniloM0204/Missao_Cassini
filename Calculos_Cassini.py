import numpy as np
import matplotlib.pyplot as plt
import os

mu_Sol = 0.0002959122 # UA^3/dia^2

mu_planetas = {
    'Venus': mu_Sol / 408523.7,
    'Terra': mu_Sol / 332946.0,
    'Jupiter': mu_Sol / 1047.3486,
    'Saturno': mu_Sol / 3498.0
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

def intervalo_SOI(dados_cassini, dados_planeta_lista, mu_p, idx_c_centro, idx_p_centro):
    """
    Encontra entrada/saída da SOI sincronizando os índices da Cassini e do Planeta.
    """
    #Cálculo do Raio da SOI
    massa_rel = mu_p / mu_Sol
    r_planeta_sol = np.linalg.norm(dados_planeta_lista[idx_p_centro]['r'])
    r_soi = r_planeta_sol * (massa_rel ** 0.4)
    offset = idx_p_centro - idx_c_centro

    n_c = len(dados_cassini)
    n_p = len(dados_planeta_lista)

    idx_pre = idx_c_centro
    idx_pos = idx_c_centro

    #Busca para Entrada
    while idx_pre > 0:
        idx_p_atual = idx_pre + offset
        if idx_p_atual < 0 or idx_p_atual >= n_p:
            break

        dist = np.linalg.norm(dados_cassini[idx_pre]['r'] - dados_planeta_lista[idx_p_atual]['r'])
        if dist > r_soi: # entrando na SOI
            break
        idx_pre -= 1

    #Busca para Saída
    while idx_pos < n_c - 1:
        idx_p_atual = idx_pos + offset
        if idx_p_atual < 0 or idx_p_atual >= n_p:
            break

        dist = np.linalg.norm(dados_cassini[idx_pos]['r'] - dados_planeta_lista[idx_p_atual]['r'])
        if dist > r_soi: # Saiu da SOI
            break
        idx_pos += 1

    if idx_pos - idx_pre < 4:
        idx_pre = max(0, idx_c_centro - 2)
        idx_pos = min(n_c - 1, idx_c_centro + 2)

    return idx_pre, idx_pos, r_soi


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

    # Raio do anel
    raio_anel = 0.25

    x_anel = xs_real + raio_anel * np.cos(theta_anel)
    y_anel = ys_real + (raio_anel) * np.sin(theta_anel)

    # Plotar o anel
    plt.plot(x_anel, y_anel, '-', color='peru', alpha=0.8, linewidth=1.5)

plt.legend(loc='upper left', fontsize='small', framealpha=0.9, bbox_to_anchor=(1, 1))
plt.tight_layout()
plt.xlim(-2,10)
plt.ylim(-4,6)
plt.savefig("./Outputs/trajetoria_completa_cassini.png", dpi=300)
plt.show()

# Calculo Flybys e Trechos
ocorre_flyby = [('Venus1', 'Venus'), ('Venus2', 'Venus'), ('Terra', 'Terra'), ('Jupiter', 'Jupiter')]

with open("./Outputs/resultados_flyby.txt", "w") as f:
    f.write("Resultados Analiticos dos Flybys pela SOI:\n\n")

for flyby, planeta in ocorre_flyby:
    jd_flyby = eventos_jd[flyby]

    idx_flyby_c = (np.abs(jds - jd_flyby)).argmin() # Índice Cassini

    jds_planeta_arr = np.array([d['jd'] for d in dados_planetas[planeta]])
    idx_flyby_p = (np.abs(jds_planeta_arr - jd_flyby)).argmin() # Índice Planeta

    idx_pre, idx_pos, raio_soi = intervalo_SOI(
        dados_cassini,
        dados_planetas[planeta],
        mu_planetas[planeta],
        idx_flyby_c,
        idx_flyby_p
    )

    #Extrai dados dos pontos de interesse
    ponto_flyby = dados_cassini[idx_flyby_c] # Periastro
    ponto_pre = dados_cassini[idx_pre]       # Entrada SOI
    ponto_pos = dados_cassini[idx_pos]       # Saida SOI

# Planetas nos momentos correspondentes
    offset = idx_flyby_p - idx_flyby_c

    # Índices alvo para o planeta
    idx_p_pre_target = idx_pre + offset
    idx_p_pos_target = idx_pos + offset

    tam_dados_planeta = len(dados_planetas[planeta])

    # Verifica o limite no pre
    if idx_p_pre_target < 0:
        idx_p_pre_target = 0
        idx_pre = idx_p_pre_target - offset # Ajusta Cassini para o mesmo dia

    # Verifica o limite no pos
    if idx_p_pos_target >= tam_dados_planeta:
        idx_p_pos_target = tam_dados_planeta - 1
        idx_pos = idx_p_pos_target - offset # Ajusta Cassini para o mesmo dia

    # Acessa os dados
    ponto_planeta_exact = dados_planetas[planeta][idx_flyby_p]
    ponto_planeta_pre   = dados_planetas[planeta][idx_p_pre_target]
    ponto_planeta_pos   = dados_planetas[planeta][idx_p_pos_target]

    # Atualiza os pontos da Cassini caso os índices tenham sido ajustados acima
    ponto_pre = dados_cassini[idx_pre]
    ponto_pos = dados_cassini[idx_pos]

    # Calculo do Flyby teorico
    orb_momento = Orbita()
    orb_momento.vet_estado(mu_Sol, ponto_flyby['r'], ponto_flyby['v'])
    r_p, v_inf, e_h, delta = orb_momento.flyby(ponto_planeta_exact['r'], ponto_planeta_exact['v'], mu_planetas[planeta])

    ganho = orb_momento.ganho_flyby(v_inf, delta)
    ganho_km = ganho * 1731.45683
    v_inf_km = v_inf * 1731.45683

    # Calculo do Flyby observado nas bordas da SOI
    r_rel_pre = ponto_pre['r'] - ponto_planeta_pre['r']

    v_rel_pre = ponto_pre['v'] - ponto_planeta_pre['v']
    v_rel_pos = ponto_pos['v'] - ponto_planeta_pos['v']

    mag_pre_km = np.linalg.norm(v_rel_pre) * 1731.45683
    mag_pos_km = np.linalg.norm(v_rel_pos) * 1731.45683

    # Excentricidade na SOI
    orb_soi = Orbita()
    orb_soi.vet_estado(mu_planetas[planeta], r_rel_pre, v_rel_pre)
    e_observado = orb_soi.e

    # Deflexão
    cos_theta = np.dot(v_rel_pre, v_rel_pos) / (np.linalg.norm(v_rel_pre) * np.linalg.norm(v_rel_pos))
    cos_theta = np.clip(cos_theta, -1.0, 1.0)
    deflexao_observada = np.degrees(np.arccos(cos_theta))


    with open("./Outputs/resultados_flyby.txt", "a") as f:
        f.write(f"FLYBY: {flyby} ({planeta})\n")
        f.write("\n")
        f.write(" Flyby Teorico\n")
        f.write(f"    - Excentricidade: {e_h:.4f}\n")
        f.write(f"    - Deflexao Max:   {delta:.2f} graus\n")
        f.write(f"    - Delta V:  {ganho_km:.2f} km/s\n")
        f.write("\n")
        f.write(" Flyby na SOI\n")
        f.write(f"    - Excentricidade SOI: {e_observado:.4f}\n")
        f.write(f"    - V_rel Entrada:  {mag_pre_km:.2f} km/s\n")
        f.write(f"    - V_rel Saida:    {mag_pos_km:.2f} km/s\n")
        f.write(f"    - Deflexao Real:  {deflexao_observada:.2f} graus\n")
        f.write("\n")
