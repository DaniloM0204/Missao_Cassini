import numpy as np
import matplotlib.pyplot as plt
import os

μ_Sol = 0.0002959122 # AU^3/dia^2

eventos_jd = {
    'Lança_Terra': 2450737.5,  # Out 1997
    'Venus1':     2450929.5,  # Abr 1998
    'Venus2':     2451353.5,  # Jun 1999
    'Terra':      2451408.5,  # Ago 1999
    'Jupiter':    2451908.5,  # Dez 2000
    'Saturno':    2453187.5   # Jul 2004
}

cores_fases = ["blue", 'orange', 'green', 'red', 'purple']

def rotacionar_90_horario(x, y):
    x_rot = y
    y_rot = -x
    return x_rot, y_rot


class Orbita:
    def __init__(self):
        self.secao_nome = ""
        self.a = 0
        self.e = 0
        self.p = 0
        self.omega = 0
        self.rvec = np.zeros(3)
        self.vvec = np.zeros(3)
        self.e_vec = np.zeros(3)

    def vet_estado(self, μ_Sol, rvec, vvec):
        self.rvec = rvec
        self.vvec = vvec

        r = np.linalg.norm(rvec)
        v = np.linalg.norm(vvec)

        self.a = 1 / ((2 / r) - (v**2 / μ_Sol))

        hvec = np.cross(rvec, vvec)

        termo1 = (v**2 - (μ_Sol / r)) * rvec
        termo2 = np.dot(rvec, vvec) * vvec
        self.e_vec = (1 / μ_Sol) * (termo1 - termo2)

        self.e = np.linalg.norm(self.e_vec)
        self.p = np.linalg.norm(hvec)**2 / μ_Sol

        self.omega = np.arctan2(self.e_vec[1], self.e_vec[0])

    def raio_orbita(self, f):
        return (self.p) / (1 + self.e * np.cos(f))

    def pos_orb(self, ω, f):
        xpq = self.raio_orbita(f) * np.cos(ω + f)
        ypq = self.raio_orbita(f) * np.sin(ω + f)
        return xpq, ypq

def ler_horizons_completo(caminho_arquivo):
    """
    Lê JD, X, Y, Z, VX, VY, VZ dos arquivos .txt
    """
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
                        # Índices típicos do Horizons (JD, Date, X, Y, Z, VX, VY, VZ)
                        jd = float(parts[0])
                        r = np.array([float(parts[2]), float(parts[3]), float(parts[4])])
                        v = np.array([float(parts[5]), float(parts[6]), float(parts[7])])
                        dados.append({'jd': jd, 'r': r, 'v': v})
                    except (ValueError, IndexError):
                        continue
    except FileNotFoundError:
        print(f"Arquivo não encontrado: {caminho_arquivo}")
        return []
    return dados

arquivos_planetas = {
    'Venus': 'VetoresVenus.txt',
    'Terra': 'VetoresTerra.txt',
    'Jupiter': 'VetoresJupiter.txt',
    'Saturno': 'VetoresSaturno.txt'
}
dados_planetas = {}
for nome, arquivo in arquivos_planetas.items():
    caminho = os.path.join("Dados", arquivo)
    dados = ler_horizons_completo(caminho)
    if dados:
        dados_planetas[nome] = dados

dados_cassini = ler_horizons_completo(os.path.join("Dados", "VetoresCassini.txt"))

plt.figure(figsize=(12, 12))
plt.title("Trajetória da Cassini")
plt.xlabel("Distância [AU]")
plt.ylabel("Distância [AU]")
plt.axis('equal')
plt.grid(True, linestyle='--', alpha=0.3)

plt.plot(0, 0, 'yo', markersize=12, label='Sol')

#Nodo Ascendente
plt.arrow(8, 0, 2, 0, color='black', alpha=0.6, width=0.02, head_width=0.1)
plt.text(10, -0.5, "Nodo γ", fontsize=10, verticalalignment='center', fontweight='bold')

f_vals = np.linspace(0, 2 * np.pi, 300)
cor_planeta = {'Venus': 'gray', 'Terra': 'blue', 'Jupiter': 'brown', 'Saturno': 'gold'}


if 'dados_planetas' in locals():
    for nome, dados in dados_planetas.items():
        ponto = dados[0]
        orb = Orbita()
        orb.vet_estado(μ_Sol, ponto['r'], ponto['v'])

        # Calcula órbita
        x_elipse, y_elipse = orb.pos_orb(orb.omega, f_vals)

        x_rot, y_rot = rotacionar_90_horario(x_elipse, y_elipse)
        # plt.plot(x_rot, y_rot, '--', color=cor_planeta.get(nome, 'gray'), alpha=0.3, linewidth=1)

        # px_rot, py_rot = rotacionar_90_horario(ponto['r'][0], ponto['r'][1])
        # plt.plot(px_rot, py_rot, 'o', color=cor_planeta.get(nome, 'gray'), markersize=4)

fases = [
    ('Terra -> Venus 1', eventos_jd['Lança_Terra'], eventos_jd['Venus1']),
    ('Venus 1 -> Venus 2',    eventos_jd['Venus1'],     eventos_jd['Venus2']),
    ('Venus 2 -> Terra',      eventos_jd['Venus2'],     eventos_jd['Terra']),
    ('Terra -> Jupiter',      eventos_jd['Terra'],      eventos_jd['Jupiter']),
    ('Jupiter -> Saturno',    eventos_jd['Jupiter'],    eventos_jd['Saturno'])
]

if dados_cassini:
    jds = np.array([d['jd'] for d in dados_cassini])

    for i, (nome_fase, jd_inicio, jd_fim) in enumerate(fases):
        cor = cores_fases[i % len(cores_fases)]

        indices = np.where((jds >= jd_inicio) & (jds <= jd_fim))[0]

        if len(indices) > 0:
            #Vetores de posição e velocidade
            segmento_r = np.array([dados_cassini[k]['r'] for k in indices])
            segmento_v = np.array([dados_cassini[k]['v'] for k in indices])

            #Elipse de transferência
            idx_mid = indices[len(indices)//4]
            ponto_ref = dados_cassini[idx_mid]

            orb_fase = Orbita()
            orb_fase.vet_estado(μ_Sol, ponto_ref['r'], ponto_ref['v'])

            #Calcular Elipse
            x_elipse, y_elipse = orb_fase.pos_orb(orb_fase.omega, f_vals)
            x_elipse_rot, y_elipse_rot = rotacionar_90_horario(x_elipse, y_elipse)

            plt.plot(x_elipse_rot, y_elipse_rot, ':', color=cor, alpha=0.5, linewidth=1)

            #Rotacionar Caminho Real
            x_real = segmento_r[:, 0]
            y_real = segmento_r[:, 1]

            x_real_rot, y_real_rot = rotacionar_90_horario(x_real, y_real)

            plt.plot(x_real_rot, y_real_rot, '-', color=cor, linewidth=1.5, label=nome_fase)

            plt.plot(x_real_rot[0], y_real_rot[0], 'o', color=cor, markersize=4)

plt.legend(loc='upper left', fontsize='small', framealpha=0.9, bbox_to_anchor=(1, 1))
plt.tight_layout()
plt.savefig("./Plot_Orbitas/trajetoria_completa_cassini.png", dpi=300)
plt.show()
