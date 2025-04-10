# Kaos
Kaos Teorisi Simülatörü
Collecting workspace information# Kaos Teorisi Simülatörü

Bu proje, kaos teorisinin temel kavramlarını interaktif bir şekilde keşfetmek için tasarlanmış bir Python uygulamasıdır. Kullanıcı arayüzü üzerinden parametre değerlerini değiştirerek kaotik sistemlerin davranışlarını gerçek zamanlı olarak gözlemleyebilirsiniz.

## Özellikler

- **Lorenz Çekicisi**: Atmosferik konveksiyonu modelleyen ünlü kaotik sistemin 3D simülasyonu
- **Çift Sarkaç**: Kaotik hareketin klasik bir örneği olan çift sarkaç sistemi ve faz uzayı gösterimi
- **Lojistik Harita**: Popülasyon dinamiklerini modelleyen basit ama zengin kaotik davranışlar sergileyen lojistik haritanın zaman serisi ve çatallanma diyagramı

## Gereksinimler

- Python 3.6+
- NumPy
- SciPy
- Matplotlib
- PyQt5

## Kurulum


# Gerekli kütüphaneleri yükleyin
pip install numpy scipy matplotlib PyQt5


## Kullanım

python analog.py


## Sistemler Hakkında

### Lorenz Çekicisi

Edward Lorenz tarafından 1963 yılında keşfedilen bu kaotik sistem, atmosferik konveksiyonu basitleştirilmiş bir şekilde modellemektedir. Üç parametreye sahiptir:

- σ (sigma): Prandtl sayısı
- ρ (rho): Rayleigh sayısı 
- β (beta): Sistemin geometrik özelliğiyle ilgili parametre

Başlangıç koşullarına hassas bağımlılık gösteren klasik bir kaotik sistemdir.

### Çift Sarkaç

Çift sarkaç, kaotik davranışı gösteren en basit mekanik sistemlerden biridir. İki sarkaç birbirine bağlanmıştır ve hareket denklemleri doğrusal olmayan terimler içerir. Çift sarkacın başlangıç koşullarındaki küçük değişiklikler, zamanla tamamen farklı davranışlara yol açar.

### Lojistik Harita

Lojistik harita, nüfus dinamiklerini modellemek için kullanılan basit bir yinelemeli fonksiyondur. `x(n+1) = r * x(n) * (1 - x(n))` denklemiyle tanımlanır. Parametre `r`'nin değerine bağlı olarak:

- r < 3.0: Sistem tek bir sabit noktaya yakınsar
- 3.0 < r < 3.57: Sistem periyodik davranış gösterir (periyot katlamaları)
- r > 3.57: Sistem kaotik davranış sergiler (bazı düzenli aralıklar hariç)

Çatallanma diyagramı, parametre `r` değiştikçe sistemin davranışının nasıl değiştiğini gösterir.

## Katkıda Bulunma

Öneriler, hata bildirimleri ve katkılar için lütfen bir Issue açın veya Pull Request gönderin.

Bu simülatör, kaos teorisinin eğitim amaçlı anlaşılmasını kolaylaştırmak için geliştirilmiştir. Sistemlerin davranışlarını keşfetmek için parametre değerlerini değiştirerek deneyler yapabilirsiniz.
