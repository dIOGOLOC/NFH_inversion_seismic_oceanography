import numpy as np
import pyvista as pv

'''
Code to generate a moving gif from an active plotter.


Inspired by https://tutorial.pyvista.org/tutorial/03_figures/d_gif.html
'''
def pyvista_create_gif(folder_input):
    block = folder_input+'/DATA/density_block.vtk'
    
    grid = pv.read(block)
    
    # Configurar o plotter
    plotter = pv.Plotter(off_screen=True)  # off_screen=True para rodar sem janela aberta
    
    # Adicionar o mesh com limites de cor e bordas
    plotter.add_mesh(
        grid,
        scalars="values",
        cmap="cividis",
        clim=(1.00, 2.00),
        show_edges=True,
        show_scalar_bar=False
    )
    
    
    # Adicionar manualmente a barra de cores com título
    plotter.add_scalar_bar(
        title="Density (g/cm³)",  # Rótulo da barra de cores
        title_font_size=12,        # Tamanho do texto do rótulo
        label_font_size=10,        # Tamanho do texto dos números
        fmt="%.0f",                # Formato dos números
        position_x=0.25,           # Ajuste horizontal (0.0 à esquerda, 1.0 à direita)
        position_y=0.9,            # Ajuste vertical (0.0 embaixo, 1.0 em cima)
        width=0.5,                 # Largura relativa da barra de cores
        height=0.08,               # Altura relativa da barra de cores
        vertical=False             # Barra de cores horizontal  
    )
    
    # Configurar GIF
    gif_filename = folder_input+"/FIGURES/block_soil_model.gif"
    plotter.open_gif(gif_filename)
    
    # Criar a rotação ao redor do eixo central
    n_frames = 60  # Número de quadros do GIF
    distance = 100  # Aumentar a distância para visualizar todo o bloco
    inclination_angle = 30  # Ângulo de inclinação em graus
    
    # Centro do bloco para a rotação
    center = (25, 25,12)
    
    for i, angle in enumerate(np.linspace(0, 360, n_frames)):
        # Calcular a posição da câmera ao redor do centro do bloco
        z_inclination = np.sin(np.radians(inclination_angle))  # Inclinação no eixo Z
        x_pos = center[0] + distance * np.cos(np.radians(angle))
        y_pos = center[1] + distance * np.sin(np.radians(angle))
        z_pos = center[2] + distance * z_inclination  # Ajuste do Z para inclinação
        
        # Definir a posição da câmera com base no ângulo e centro
        plotter.camera_position = [
            (x_pos, y_pos, z_pos),  # Posição da câmera ao redor do centro
            center,                 # Foco da câmera (centro do bloco)
            (0, 0, 1),              # Vetor "up"
        ]
        
        plotter.write_frame()  # Escreve o quadro no GIF
    
    plotter.close()