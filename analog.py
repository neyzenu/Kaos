import sys
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use('Qt5Agg')
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from PyQt5.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, 
                             QHBoxLayout, QComboBox, QSlider, QLabel, QPushButton,
                             QGroupBox, QGridLayout, QTabWidget, QScrollArea)
from PyQt5.QtCore import Qt
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

class KaosSimulatoru(QMainWindow):
    def __init__(self):
        super().__init__()
        
        # Ana pencere ayarları
        self.setWindowTitle("Kaos Teorisi Simülatörü")
        self.setGeometry(100, 100, 1200, 800)
        self.setMinimumSize(800, 600)  # Minimum pencere boyutu ekle
        
        # Ana widget ve layout
        self.central_widget = QWidget()
        self.setCentralWidget(self.central_widget)
        self.layout = QHBoxLayout(self.central_widget)
        
        # Tab widget oluşturma
        self.tab_widget = QTabWidget()
        
        # İlk sekme değişikliğini önle (grafik alanları oluşturulmadan önce)
        self.ignore_first_tab_change = True
        self.tab_widget.currentChanged.connect(self.tab_changed)
        
        # Sistem seçimi için kaydırılabilir kontrol paneli
        self.kontrol_scroll = QScrollArea()
        self.kontrol_scroll.setWidgetResizable(True)
        self.kontrol_scroll.setMinimumWidth(300)  # Kontrol paneli minimum genişliği
        self.kontrol_scroll.setMaximumWidth(400)  # Kontrol paneli maximum genişliği
        
        self.kontrol_panel = QWidget()
        self.kontrol_layout = QVBoxLayout(self.kontrol_panel)
        self.kontrol_layout.setContentsMargins(8, 8, 8, 8)
        self.kontrol_layout.setSpacing(10)
        
        self.kontrol_scroll.setWidget(self.kontrol_panel)
        
        # Sekmeler
        self.lorenz_tab = QWidget()
        self.double_pendulum_tab = QWidget()
        self.logistic_tab = QWidget()
        
        # Sekmeleri hazırla
        self.setup_lorenz_tab()
        self.setup_double_pendulum_tab()
        self.setup_logistic_tab()
        
        self.tab_widget.addTab(self.lorenz_tab, "Lorenz Çekicisi")
        self.tab_widget.addTab(self.double_pendulum_tab, "Çift Sarkaç")
        self.tab_widget.addTab(self.logistic_tab, "Lojistik Harita")
        
        # Ana düzeni tamamla
        self.layout.addWidget(self.kontrol_scroll, 1)
        self.layout.addWidget(self.tab_widget, 3)
        
        # İlk simülasyonu başlat
        self.current_system = "lorenz"
        self.setup_lorenz_controls()
        self.run_simulation()

    def setup_lorenz_tab(self):
        layout = QVBoxLayout(self.lorenz_tab)
        layout.setContentsMargins(0, 0, 0, 0)
        
        # Lorenz grafik
        self.lorenz_fig = Figure(figsize=(10, 8))
        self.lorenz_canvas = FigureCanvas(self.lorenz_fig)
        self.lorenz_ax = self.lorenz_fig.add_subplot(111, projection='3d')
        self.lorenz_toolbar = NavigationToolbar(self.lorenz_canvas, self)
        
        layout.addWidget(self.lorenz_toolbar)
        layout.addWidget(self.lorenz_canvas)

    def setup_double_pendulum_tab(self):
        layout = QVBoxLayout(self.double_pendulum_tab)
        layout.setContentsMargins(0, 0, 0, 0)
        
        # Çift sarkaç grafikleri
        self.pendulum_fig = Figure(figsize=(10, 8))
        self.pendulum_canvas = FigureCanvas(self.pendulum_fig)
        
        # İki grafik oluştur: biri animasyon için, diğeri faz uzayı için
        self.pendulum_ax1 = self.pendulum_fig.add_subplot(121)
        self.pendulum_ax2 = self.pendulum_fig.add_subplot(122)
        
        self.pendulum_toolbar = NavigationToolbar(self.pendulum_canvas, self)
        
        layout.addWidget(self.pendulum_toolbar)
        layout.addWidget(self.pendulum_canvas)

    def setup_logistic_tab(self):
        layout = QVBoxLayout(self.logistic_tab)
        layout.setContentsMargins(0, 0, 0, 0)
        
        # Lojistik harita grafikleri
        self.logistic_fig = Figure(figsize=(10, 8))
        self.logistic_canvas = FigureCanvas(self.logistic_fig)
        
        # İki grafik oluştur: biri zaman serisi, diğeri çatallanma diyagramı için
        self.logistic_ax1 = self.logistic_fig.add_subplot(211)  # Zaman serisi
        self.logistic_ax2 = self.logistic_fig.add_subplot(212)  # Çatallanma diyagramı
        
        self.logistic_toolbar = NavigationToolbar(self.logistic_canvas, self)
        
        layout.addWidget(self.logistic_toolbar)
        layout.addWidget(self.logistic_canvas)

    def setup_slider(self, label_text, min_val, max_val, default_val, divisor=10.0, system_type="lorenz"):
        """Standart slider ayarı için yardımcı fonksiyon"""
        label = QLabel(f"{label_text}: {default_val/divisor}")
        slider = QSlider(Qt.Horizontal)
        slider.setRange(min_val, max_val)
        slider.setValue(default_val)
        slider.valueChanged.connect(lambda: self.slider_changed(system_type))
        
        return slider, label
        
    def setup_lorenz_controls(self):
        # Lorenz parametreleri için kontrol grubu
        lorenz_group = QGroupBox("Lorenz Parametreleri")
        grid_layout = QGridLayout()
        grid_layout.setColumnStretch(1, 1)  # Slider sütunu genişletilsin
        
        # Sigma parametresi (σ)
        self.sigma_slider, self.sigma_label = self.setup_slider("σ (sigma)", 1, 300, 100)
        grid_layout.addWidget(QLabel("σ:"), 0, 0)
        grid_layout.addWidget(self.sigma_slider, 0, 1)
        grid_layout.addWidget(self.sigma_label, 0, 2)
        
        # Rho parametresi (ρ)
        self.rho_slider, self.rho_label = self.setup_slider("ρ (rho)", 1, 1000, 280)
        grid_layout.addWidget(QLabel("ρ:"), 1, 0)
        grid_layout.addWidget(self.rho_slider, 1, 1)
        grid_layout.addWidget(self.rho_label, 1, 2)
        
        # Beta parametresi (β)
        self.beta_slider, self.beta_label = self.setup_slider("β (beta)", 1, 300, 83, 30.0)
        grid_layout.addWidget(QLabel("β:"), 2, 0)
        grid_layout.addWidget(self.beta_slider, 2, 1)
        grid_layout.addWidget(self.beta_label, 2, 2)
        
        # Başlangıç koşulları
        init_group = QGroupBox("Başlangıç Koşulları")
        init_layout = QGridLayout()
        init_layout.setColumnStretch(1, 1)  # Slider sütunu genişletilsin
        
        self.x0_slider, self.x0_label = self.setup_slider("x₀", -200, 200, 10)
        init_layout.addWidget(QLabel("x₀:"), 0, 0)
        init_layout.addWidget(self.x0_slider, 0, 1)
        init_layout.addWidget(self.x0_label, 0, 2)
        
        self.y0_slider, self.y0_label = self.setup_slider("y₀", -200, 200, 10)
        init_layout.addWidget(QLabel("y₀:"), 1, 0)
        init_layout.addWidget(self.y0_slider, 1, 1)
        init_layout.addWidget(self.y0_label, 1, 2)
        
        self.z0_slider, self.z0_label = self.setup_slider("z₀", -200, 200, 10)
        init_layout.addWidget(QLabel("z₀:"), 2, 0)
        init_layout.addWidget(self.z0_slider, 2, 1)
        init_layout.addWidget(self.z0_label, 2, 2)
        
        init_group.setLayout(init_layout)
        
        # Simülasyon süresi
        time_group = QGroupBox("Simülasyon Ayarları")
        time_layout = QGridLayout()
        time_layout.setColumnStretch(1, 1)  # Slider sütunu genişletilsin
        
        self.time_slider, self.time_label = self.setup_slider("Süre", 10, 1000, 500)
        time_layout.addWidget(QLabel("Süre:"), 0, 0)
        time_layout.addWidget(self.time_slider, 0, 1)
        time_layout.addWidget(self.time_label, 0, 2)
        
        self.points_slider, self.points_label = self.setup_slider("Nokta Sayısı", 100, 20000, 5000, 1.0)
        time_layout.addWidget(QLabel("Noktalar:"), 1, 0)
        time_layout.addWidget(self.points_slider, 1, 1)
        time_layout.addWidget(self.points_label, 1, 2)
        
        time_group.setLayout(time_layout)
        
        # Tüm kontrol panelini temizle ve yeni kontrolleri ekle
        self.clear_layout_except_current_tab_controls()
        
        # Grup kutularını ana yerleşime ekle
        lorenz_group.setLayout(grid_layout)
        self.kontrol_layout.addWidget(lorenz_group)
        self.kontrol_layout.addWidget(init_group)
        self.kontrol_layout.addWidget(time_group)
        
        # Bu simülasyon çalıştırma butonu
        run_button = QPushButton("Simülasyonu Yenile")
        run_button.clicked.connect(self.run_simulation)
        self.kontrol_layout.addWidget(run_button)
        
        # Bilgi etiketi
        info_label = QLabel("Değerleri değiştirdikçe simülasyon otomatik olarak güncellenir.")
        info_label.setWordWrap(True)
        info_label.setStyleSheet("color: #666; font-style: italic;")
        self.kontrol_layout.addWidget(info_label)
        
        self.kontrol_layout.addStretch()  # Boş alanı aşağıya doğru genişlet

    def setup_pendulum_controls(self):
        # Çift sarkaç parametreleri
        pendulum_group = QGroupBox("Çift Sarkaç Parametreleri")
        grid_layout = QGridLayout()
        grid_layout.setColumnStretch(1, 1)
        
        # Sarkaç 1 kütlesi
        self.m1_slider, self.m1_label = self.setup_slider("m₁ (kütle 1)", 1, 100, 10, system_type="pendulum")
        grid_layout.addWidget(QLabel("m₁:"), 0, 0)
        grid_layout.addWidget(self.m1_slider, 0, 1)
        grid_layout.addWidget(self.m1_label, 0, 2)
        
        # Sarkaç 2 kütlesi
        self.m2_slider, self.m2_label = self.setup_slider("m₂ (kütle 2)", 1, 100, 10, system_type="pendulum")
        grid_layout.addWidget(QLabel("m₂:"), 1, 0)
        grid_layout.addWidget(self.m2_slider, 1, 1)
        grid_layout.addWidget(self.m2_label, 1, 2)
        
        # Sarkaç 1 uzunluğu
        self.l1_slider, self.l1_label = self.setup_slider("L₁ (uzunluk 1)", 1, 100, 10, system_type="pendulum")
        grid_layout.addWidget(QLabel("L₁:"), 2, 0)
        grid_layout.addWidget(self.l1_slider, 2, 1)
        grid_layout.addWidget(self.l1_label, 2, 2)
        
        # Sarkaç 2 uzunluğu
        self.l2_slider, self.l2_label = self.setup_slider("L₂ (uzunluk 2)", 1, 100, 10, system_type="pendulum")
        grid_layout.addWidget(QLabel("L₂:"), 3, 0)
        grid_layout.addWidget(self.l2_slider, 3, 1)
        grid_layout.addWidget(self.l2_label, 3, 2)
        
        # Başlangıç koşulları
        init_group = QGroupBox("Başlangıç Koşulları")
        init_layout = QGridLayout()
        init_layout.setColumnStretch(1, 1)
        
        # Başlangıç açısı 1
        self.theta1_slider, self.theta1_label = self.setup_slider("θ₁", 0, 360, 120, 1.0, "pendulum")
        init_layout.addWidget(QLabel("θ₁:"), 0, 0)
        init_layout.addWidget(self.theta1_slider, 0, 1)
        init_layout.addWidget(self.theta1_label, 0, 2)
        
        # Başlangıç açısı 2
        self.theta2_slider, self.theta2_label = self.setup_slider("θ₂", 0, 360, 180, 1.0, "pendulum")
        init_layout.addWidget(QLabel("θ₂:"), 1, 0)
        init_layout.addWidget(self.theta2_slider, 1, 1)
        init_layout.addWidget(self.theta2_label, 1, 2)
        
        # Yerçekimi
        self.g_slider, self.g_label = self.setup_slider("g (yerçekimi)", 1, 300, 98, system_type="pendulum")
        init_layout.addWidget(QLabel("g:"), 2, 0)
        init_layout.addWidget(self.g_slider, 2, 1)
        init_layout.addWidget(self.g_label, 2, 2)
        
        init_group.setLayout(init_layout)
        
        # Simülasyon ayarları
        time_group = QGroupBox("Simülasyon Ayarları")
        time_layout = QGridLayout()
        time_layout.setColumnStretch(1, 1)
        
        self.pendulum_time_slider, self.pendulum_time_label = self.setup_slider("Süre", 10, 1000, 300, system_type="pendulum")
        time_layout.addWidget(QLabel("Süre:"), 0, 0)
        time_layout.addWidget(self.pendulum_time_slider, 0, 1)
        time_layout.addWidget(self.pendulum_time_label, 0, 2)
        
        time_group.setLayout(time_layout)
        
        # Tüm kontrol panelini temizle ve yeni kontrolleri ekle
        self.clear_layout_except_current_tab_controls()
        
        # Grup kutularını ekle
        pendulum_group.setLayout(grid_layout)
        self.kontrol_layout.addWidget(pendulum_group)
        self.kontrol_layout.addWidget(init_group)
        self.kontrol_layout.addWidget(time_group)
        
        # Bu simülasyon çalıştırma butonu
        run_button = QPushButton("Simülasyonu Yenile")
        run_button.clicked.connect(self.run_simulation)
        self.kontrol_layout.addWidget(run_button)
        
        # Bilgi etiketi
        info_label = QLabel("Değerleri değiştirdikçe simülasyon otomatik olarak güncellenir. İki sarkacın farklı başlangıç koşullarındaki kaotik hareketini gözlemleyin.")
        info_label.setWordWrap(True)
        info_label.setStyleSheet("color: #666; font-style: italic;")
        self.kontrol_layout.addWidget(info_label)
        
        self.kontrol_layout.addStretch()

    def setup_logistic_controls(self):
        # Lojistik harita parametreleri
        logistic_group = QGroupBox("Lojistik Harita Parametreleri")
        grid_layout = QGridLayout()
        grid_layout.setColumnStretch(1, 1)
        
        # r parametresi
        self.r_slider, self.r_label = self.setup_slider("r", 1, 40, 35, system_type="logistic")
        grid_layout.addWidget(QLabel("r:"), 0, 0)
        grid_layout.addWidget(self.r_slider, 0, 1)
        grid_layout.addWidget(self.r_label, 0, 2)
        
        # Başlangıç koşulu
        self.x0_logistic_slider, self.x0_logistic_label = self.setup_slider("x₀", 1, 99, 50, 100.0, "logistic")
        grid_layout.addWidget(QLabel("x₀:"), 1, 0)
        grid_layout.addWidget(self.x0_logistic_slider, 1, 1)
        grid_layout.addWidget(self.x0_logistic_label, 1, 2)
        
        # Simülasyon ayarları
        self.iterations_slider, self.iterations_label = self.setup_slider("İterasyon", 10, 1000, 100, 1.0, "logistic")
        grid_layout.addWidget(QLabel("İterasyon:"), 2, 0)
        grid_layout.addWidget(self.iterations_slider, 2, 1)
        grid_layout.addWidget(self.iterations_label, 2, 2)
        
        # Çatallanma diyagramı için ayarlar
        bifurcation_group = QGroupBox("Çatallanma Diyagramı")
        bifurcation_layout = QGridLayout()
        bifurcation_layout.setColumnStretch(1, 1)
        
        self.r_min_slider, self.r_min_label = self.setup_slider("r min", 1, 39, 25, system_type="logistic")
        bifurcation_layout.addWidget(QLabel("r min:"), 0, 0)
        bifurcation_layout.addWidget(self.r_min_slider, 0, 1)
        bifurcation_layout.addWidget(self.r_min_label, 0, 2)
        
        self.r_max_slider, self.r_max_label = self.setup_slider("r max", 25, 40, 40, system_type="logistic")
        bifurcation_layout.addWidget(QLabel("r max:"), 1, 0)
        bifurcation_layout.addWidget(self.r_max_slider, 1, 1)
        bifurcation_layout.addWidget(self.r_max_label, 1, 2)
        
        bifurcation_group.setLayout(bifurcation_layout)
        
        # Tüm kontrol panelini temizle ve yeni kontrolleri ekle
        self.clear_layout_except_current_tab_controls()
        
        # Grup kutularını ekle
        logistic_group.setLayout(grid_layout)
        self.kontrol_layout.addWidget(logistic_group)
        self.kontrol_layout.addWidget(bifurcation_group)
        
        # Bu simülasyon çalıştırma butonu
        run_button = QPushButton("Simülasyonu Yenile")
        run_button.clicked.connect(self.run_simulation)
        self.kontrol_layout.addWidget(run_button)
        
        # Bilgi etiketi
        info_label = QLabel("Çatallanma diyagramı r ≈ 3.57'de başlayan periyot katlamalarını gösterir. r > 3.9 için tam kaotik rejim gözlenir.")
        info_label.setWordWrap(True)
        info_label.setStyleSheet("color: #666; font-style: italic;")
        self.kontrol_layout.addWidget(info_label)
        
        self.kontrol_layout.addStretch()

    def tab_changed(self, index):
        # İlk tab değişimini atla
        if self.ignore_first_tab_change:
            self.ignore_first_tab_change = False
            return
            
        # Sekme değiştirildiğinde çağrılır
        if index == 0:
            self.current_system = "lorenz"
            self.setup_lorenz_controls()
        elif index == 1:
            self.current_system = "pendulum"
            self.setup_pendulum_controls()
        elif index == 2:
            self.current_system = "logistic"
            self.setup_logistic_controls()
        
        self.run_simulation()

    def clear_layout_except_current_tab_controls(self):
        # Kontrol panelindeki tüm widget'ları temizle
        for i in reversed(range(self.kontrol_layout.count())): 
            widget = self.kontrol_layout.itemAt(i).widget()
            if widget is not None:
                widget.deleteLater()

    def slider_changed(self, system_type):
        # Slider değerleri değiştiğinde etiketleri güncelle ve simülasyonu çalıştır
        if system_type == "lorenz":
            self.sigma_label.setText(f"σ (sigma): {self.sigma_slider.value() / 10.0:.2f}")
            self.rho_label.setText(f"ρ (rho): {self.rho_slider.value() / 10.0:.2f}")
            self.beta_label.setText(f"β (beta): {self.beta_slider.value() / 30.0:.2f}")
            
            self.x0_label.setText(f"x₀: {self.x0_slider.value() / 10.0:.2f}")
            self.y0_label.setText(f"y₀: {self.y0_slider.value() / 10.0:.2f}")
            self.z0_label.setText(f"z₀: {self.z0_slider.value() / 10.0:.2f}")
            
            self.time_label.setText(f"Süre: {self.time_slider.value() / 10.0:.2f}")
            self.points_label.setText(f"Nokta Sayısı: {self.points_slider.value()}")
        
        elif system_type == "pendulum":
            self.m1_label.setText(f"m₁ (kütle 1): {self.m1_slider.value() / 10.0:.2f}")
            self.m2_label.setText(f"m₂ (kütle 2): {self.m2_slider.value() / 10.0:.2f}")
            self.l1_label.setText(f"L₁ (uzunluk 1): {self.l1_slider.value() / 10.0:.2f}")
            self.l2_label.setText(f"L₂ (uzunluk 2): {self.l2_slider.value() / 10.0:.2f}")
            
            self.theta1_label.setText(f"θ₁: {self.theta1_slider.value()}°")
            self.theta2_label.setText(f"θ₂: {self.theta2_slider.value()}°")
            self.g_label.setText(f"g (yerçekimi): {self.g_slider.value() / 10.0:.2f}")
            
            self.pendulum_time_label.setText(f"Süre: {self.pendulum_time_slider.value() / 10.0:.2f}")
        
        elif system_type == "logistic":
            self.r_label.setText(f"r: {self.r_slider.value() / 10.0:.2f}")
            self.x0_logistic_label.setText(f"x₀: {self.x0_logistic_slider.value() / 100.0:.2f}")
            self.iterations_label.setText(f"İterasyon: {self.iterations_slider.value()}")
            
            self.r_min_label.setText(f"r min: {self.r_min_slider.value() / 10.0:.2f}")
            self.r_max_label.setText(f"r max: {self.r_max_slider.value() / 10.0:.2f}")
        
        # Her değişiklikte simülasyonu güncelle
        self.run_simulation()

    def lorenz_system(self, t, xyz, sigma, rho, beta):
        x, y, z = xyz
        dx_dt = sigma * (y - x)
        dy_dt = x * (rho - z) - y
        dz_dt = x * y - beta * z
        return [dx_dt, dy_dt, dz_dt]

    def double_pendulum_system(self, t, y, m1, m2, l1, l2, g):
        theta1, z1, theta2, z2 = y
        
        c = np.cos(theta1 - theta2)
        s = np.sin(theta1 - theta2)
        
        # İkinci dereceden diferansiyel denklemleri çöz
        theta1_dot = z1
        z1_dot = (m2 * g * np.sin(theta2) * c - m2 * s * (l1 * z1**2 * c + l2 * z2**2) -
                 (m1 + m2) * g * np.sin(theta1)) / (l1 * (m1 + m2 * s**2))
        
        theta2_dot = z2
        z2_dot = ((m1 + m2) * (l1 * z1**2 * s - g * np.sin(theta2) + g * np.sin(theta1) * c) +
                 m2 * l2 * z2**2 * s * c) / (l2 * (m1 + m2 * s**2))
        
        return [theta1_dot, z1_dot, theta2_dot, z2_dot]

    def logistic_map(self, x0, r, iterations):
        # Lojistik harita x_(n+1) = r * x_n * (1 - x_n)
        x = np.zeros(iterations)
        x[0] = x0
        
        for i in range(iterations-1):
            x[i+1] = r * x[i] * (1 - x[i])
            
        return x

    def run_simulation(self):
        if self.current_system == "lorenz":
            self.run_lorenz_simulation()
        elif self.current_system == "pendulum":
            self.run_pendulum_simulation()
        elif self.current_system == "logistic":
            self.run_logistic_simulation()

    def run_lorenz_simulation(self):
        # Parametreleri slider değerlerinden al
        sigma = self.sigma_slider.value() / 10.0
        rho = self.rho_slider.value() / 10.0
        beta = self.beta_slider.value() / 30.0
        
        x0 = self.x0_slider.value() / 10.0
        y0 = self.y0_slider.value() / 10.0
        z0 = self.z0_slider.value() / 10.0
        
        t_max = self.time_slider.value() / 10.0
        points = self.points_slider.value()
        
        # Zaman dizisi
        t = np.linspace(0, t_max, points)
        
        # Lorenz sistemini çöz
        sol = solve_ivp(
            lambda t, xyz: self.lorenz_system(t, xyz, sigma, rho, beta),
            [0, t_max],
            [x0, y0, z0],
            t_eval=t
        )
        
        # Grafiği temizle ve çiz
        self.lorenz_ax.clear()
        self.lorenz_ax.plot(sol.y[0], sol.y[1], sol.y[2], lw=0.5)
        self.lorenz_ax.set_xlabel('X')
        self.lorenz_ax.set_ylabel('Y')
        self.lorenz_ax.set_zlabel('Z')
        self.lorenz_ax.set_title(f"Lorenz Çekicisi (σ={sigma}, ρ={rho}, β={beta:.2f})")
        
        # Grafiği güncelle
        self.lorenz_canvas.draw()

    def run_pendulum_simulation(self):
        # Parametreleri slider değerlerinden al
        m1 = self.m1_slider.value() / 10.0
        m2 = self.m2_slider.value() / 10.0
        l1 = self.l1_slider.value() / 10.0
        l2 = self.l2_slider.value() / 10.0
        
        theta1 = np.radians(self.theta1_slider.value())
        theta2 = np.radians(self.theta2_slider.value())
        g = self.g_slider.value() / 10.0
        
        t_max = self.pendulum_time_slider.value() / 10.0
        
        # Başlangıç koşulları [theta1, theta1_dot, theta2, theta2_dot]
        y0 = [theta1, 0, theta2, 0]
        
        # Zaman dizisi
        t = np.linspace(0, t_max, 1000)
        
        # Diferansiyel denklem sistemini çöz
        sol = solve_ivp(
            lambda t, y: self.double_pendulum_system(t, y, m1, m2, l1, l2, g),
            [0, t_max],
            y0,
            t_eval=t,
            method='RK45',
            rtol=1e-8
        )
        
        # Sarkacın kartezyen koordinatları
        theta1 = sol.y[0]
        theta2 = sol.y[2]
        
        x1 = l1 * np.sin(theta1)
        y1 = -l1 * np.cos(theta1)
        
        x2 = x1 + l2 * np.sin(theta2)
        y2 = y1 - l2 * np.cos(theta2)
        
        # İlk grafik - animasyon son karesini göster
        self.pendulum_ax1.clear()
        self.pendulum_ax1.plot([0, x1[-1]], [0, y1[-1]], 'k-', lw=2)
        self.pendulum_ax1.plot([x1[-1], x2[-1]], [y1[-1], y2[-1]], 'k-', lw=2)
        self.pendulum_ax1.plot(x1[-1], y1[-1], 'bo', markersize=10)
        self.pendulum_ax1.plot(x2[-1], y2[-1], 'ro', markersize=10)
        self.pendulum_ax1.plot(x1, y1, 'b-', alpha=0.3)  # İlk kütle yörüngesi
        self.pendulum_ax1.plot(x2, y2, 'r-', alpha=0.3)  # İkinci kütle yörüngesi
        
        # Grafik ayarları
        self.pendulum_ax1.set_xlim(-(l1+l2+0.5), (l1+l2+0.5))
        self.pendulum_ax1.set_ylim(-(l1+l2+0.5), (l1+l2+0.5))
        self.pendulum_ax1.set_aspect('equal')
        self.pendulum_ax1.grid(True)
        self.pendulum_ax1.set_title('Çift Sarkaç Hareketi')
        
        # İkinci grafik - faz uzayı (theta1 vs theta2)
        self.pendulum_ax2.clear()
        self.pendulum_ax2.plot(theta1, theta2, 'g-', alpha=0.8)
        self.pendulum_ax2.set_xlabel('θ₁ (rad)')
        self.pendulum_ax2.set_ylabel('θ₂ (rad)')
        self.pendulum_ax2.set_title('Faz Uzayı (θ₁ vs θ₂)')
        self.pendulum_ax2.grid(True)
        
        # Grafik güncelle
        self.pendulum_fig.tight_layout()
        self.pendulum_canvas.draw()

    def run_logistic_simulation(self):
        # Parametreleri al
        r = self.r_slider.value() / 10.0
        x0 = self.x0_logistic_slider.value() / 100.0
        iterations = self.iterations_slider.value()
        
        # Lojistik haritayı hesapla
        x = self.logistic_map(x0, r, iterations)
        
        # İlk grafik - zaman serisi
        self.logistic_ax1.clear()
        self.logistic_ax1.plot(range(iterations), x, 'b-', lw=1)
        self.logistic_ax1.set_xlabel('İterasyon')
        self.logistic_ax1.set_ylabel('x değeri')
        self.logistic_ax1.set_title(f'Lojistik Harita Zaman Serisi (r = {r})')
        self.logistic_ax1.grid(True)
        
        # İkinci grafik - çatallanma diyagramı
        r_min = self.r_min_slider.value() / 10.0
        r_max = self.r_max_slider.value() / 10.0
        
        r_range = np.linspace(r_min, r_max, 1000)
        
        # Her r değeri için son 50 iterasyonu sakla
        bifurcation_x = []
        bifurcation_r = []
        
        transient = 100  # Geçici davranışları atlamak için iterasyon sayısı
        
        for r_value in r_range:
            x_current = x0
            # Geçici davranışları atla
            for _ in range(transient):
                x_current = r_value * x_current * (1 - x_current)
            
            # Son değerleri kaydet
            for _ in range(50):
                x_current = r_value * x_current * (1 - x_current)
                bifurcation_x.append(x_current)
                bifurcation_r.append(r_value)
        
        self.logistic_ax2.clear()
        self.logistic_ax2.plot(bifurcation_r, bifurcation_x, 'k.', markersize=0.5)
        self.logistic_ax2.set_xlabel('r parametresi')
        self.logistic_ax2.set_ylabel('x değerleri')
        self.logistic_ax2.set_title('Çatallanma Diyagramı')
        
        # Grafiği güncelle
        self.logistic_fig.tight_layout()
        self.logistic_canvas.draw()

if __name__ == "__main__":
    # Yüksek DPI ekranlar için ölçekleme sorunlarını çöz 
    QApplication.setAttribute(Qt.AA_EnableHighDpiScaling, True)
    QApplication.setAttribute(Qt.AA_UseHighDpiPixmaps, True)
    
    app = QApplication(sys.argv)
    window = KaosSimulatoru()
    window.show()
    sys.exit(app.exec_())