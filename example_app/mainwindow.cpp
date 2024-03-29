#include "mainwindow.h"
#include "./ui_mainwindow.h"

#include <QChart>
#include <QChartView>
#include <QLineSeries>
#include <QLogValueAxis>
#include <QValueAxis>
#include <QtMath>

#include <atmosGOST_R_25645_166_2004.h>
extern "C" {
    #include <nrlmsise-00.h>
}

void MainWindow::SetupChart()
{
    m_chart1 = new QChart();

    m_chart1->setTitle("ГОСТ Р 25645.166-2004");
    m_chart1->legend()->setVisible(true);
    m_chart1->legend()->setAlignment(Qt::AlignRight);

    m_chart1->setTitleFont(QFont("Arial", 12, QFont::Bold));
    m_chart1->setFont(QFont("Arial", 10, QFont::Normal));


    QPen gridPen;
    gridPen.setColor("gray");
    gridPen.setWidth(1);
    gridPen.setStyle(Qt::SolidLine);

    QPen minorGridPen;
    minorGridPen.setColor("lightGray");
    minorGridPen.setWidth(1);
    minorGridPen.setStyle(Qt::DotLine);

    m_axisX = new QValueAxis();
    m_axisX->setTitleText("Высота, км");
    m_axisX->setLabelFormat("%0.0f");
    m_axisX->setGridLineVisible(true);
    m_axisX->setTickType(QValueAxis::TicksDynamic);
    m_axisX->setTickAnchor(100);
    m_axisX->setTickInterval(100);
    m_axisX->setMinorTickCount(4);
    m_axisX->setGridLineVisible(true);
    m_axisX->setGridLinePen(gridPen);
    m_axisX->setMinorGridLineVisible(true);
    m_axisX->setMinorGridLinePen(minorGridPen);
    m_axisX->setTitleFont(QFont("Arial", 10, QFont::Bold));
    m_axisX->setLabelsFont(QFont("Arial", 10, QFont::Normal));

    m_axisY = new QLogValueAxis();
    m_axisY->setTitleText("Плотность, кг/м<sup>3<sup/>");
    m_axisY->setLabelFormat("%0.2e");
    m_axisY->setBase(10);
    m_axisY->setGridLineVisible(true);
    m_axisY->setMinorTickCount(5);
    m_axisY->setGridLinePen(gridPen);
    m_axisY->setMinorGridLineVisible(true);
    m_axisY->setMinorGridLinePen(minorGridPen);
    m_axisY->setTitleFont(QFont("Arial", 10, QFont::Bold));
    m_axisY->setLabelsFont(QFont("Arial", 10, QFont::Normal));

    m_chart1->addAxis(m_axisX, Qt::AlignBottom);
    m_chart1->addAxis(m_axisY, Qt::AlignLeft);
    ui->chartView1->setChart(m_chart1);
}

void MainWindow::SetupGUI()
{
    QObject::connect(ui->chkBunchOfCurves, &QCheckBox::toggled,
                     [=](bool state){
                         ui->gbManualParameters->blockSignals(true);
                         ui->gbManualParameters->setChecked(!state);
                         ui->gbManualParameters->blockSignals(false);
                         if(state){
                             PlotBunchOfCurves();
                         }
                     });

    QObject::connect(ui->gbManualParameters, &QGroupBox::toggled,
                     [=](bool state){
                         ui->chkBunchOfCurves->blockSignals(true);
                         ui->chkBunchOfCurves->setChecked(!state);
                         ui->chkBunchOfCurves->blockSignals(false);
                         if(state){
                             PlotByParameters();
                         }
                     });

    QList<QDoubleSpinBox*> controlsToBind =
        ui->gbCurve1->findChildren<QDoubleSpinBox*>()
        + ui->gbCurve2->findChildren<QDoubleSpinBox*>();
    for(const QDoubleSpinBox * item : controlsToBind)
        QObject::connect(item, &QDoubleSpinBox::valueChanged,
                         this, &MainWindow::PlotByParameters);

    QObject::connect(ui->gbCurve1, &QGroupBox::toggled,
                     this, &MainWindow::PlotByParameters);
    QObject::connect(ui->gbCurve2, &QGroupBox::toggled,
                     this, &MainWindow::PlotByParameters);
    QObject::connect(ui->chkNrlmsise00on, &QCheckBox::toggled,
                     this, &MainWindow::PlotByParameters);
    QObject::connect(ui->chkGOSTon, &QCheckBox::toggled,
                     this, &MainWindow::PlotByParameters);
}

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    SetupChart();
    SetupGUI();
    PlotByParameters();
}

MainWindow::~MainWindow()
{
    delete ui;
}


/**
 * @brief Выполняет расчёт семейства кривых, варьируя F107 и высоту,
 * при этом F81 принимается равным F107, геоцентрическая координата X
 * изменяется с высотой а координаты Y, Z и прочие параметры принимаютсяэ
 * равными нулю
 */
void MainWindow::PlotBunchOfCurves()
{
    m_chart1->removeAllSeries();

    double F107 = 75;
    double F107_step = 25;

    while(F107 <= 250){
        double h_km = 120;
        double h_step = 20;
        double F81 = F107;
        QLineSeries * series = new QLineSeries();
        series->setName(QString("F_107 = %2")
                            .arg(F107, 0, 'f', 0));

        while(h_km <= 1500){
            double X[3] = {6371.0+h_km, 0.0, 0.0};
            double rho = atmosGOST_R_25645_166_2004(
                h_km, F107, 0.0, F81,1, X, 0.0, 0.0, 0.0, 0.0);
            series->append(h_km, rho);
            h_km += h_step;
        }
        m_chart1->addSeries(series);
        qDebug() << series->attachAxis(m_axisY);
        qDebug() << series->attachAxis(m_axisX);

        F107 += F107_step;
    }
}


/**
 * @brief Функция для преобразования геоцентрических координат XYZ
 * в геодезических координаты: широтыу (B), долготыу (L) и высоту (H)
 * Используются формулы П1.4-П1.19 из ПЗ-90.11
 */
void BLHtoXYZ(
    const double B_rad,
    const double L_rad,
    const double H_m,
    double & X_m,
    double & Y_m,
    double & Z_m)
{
    double a_m = 6378136; // большая полуось эллипсоида, ПЗ-90.11 Таблица 2
    double alpha = 1.0/298.25784;// Сжатие общеземного эллипсоида, ПЗ-90.11 Таблица 2
    // double e2 = 0.0066943662; // квадрат эксцентриситет эллипсоида, ПЗ-90.11, Таблица 3
    double e2 = (2 - alpha)*alpha; // квадрат эксцентриситет эллипсоида, П1.3
    double N = a_m * std::sqrt(1 - e2 * std::pow(std::sin(B_rad), 2) ); // радиус кривизны первого вертикала, ПЗ-90.11 П1.2
    X_m = (N+H_m)*std::cos(B_rad)*std::cos(L_rad);
    Y_m = (N+H_m)*std::cos(B_rad)*std::sin(L_rad);
    Z_m = ((1-e2)*N + H_m)*std::sin(B_rad);
    return;
}


/**
 * @brief Векторизованный вариант функции
 */
void BLHtoXYZ(
    const double BLH[3],
    double XYZ_m[3])
{
    double a_m = 6378136; // большая полуось эллипсоида, ПЗ-90.11 Таблица 2
    double alpha = 1.0/298.25784;// Сжатие общеземного эллипсоида, ПЗ-90.11 Таблица 2
    // double e2 = 0.0066943662; // квадрат эксцентриситет эллипсоида, ПЗ-90.11, Таблица 3
    double e2 = (2 - alpha)*alpha; // квадрат эксцентриситет эллипсоида, П1.3
    double N = a_m * std::sqrt(1 - e2 * std::pow(std::sin(BLH[0]), 2) ); // радиус кривизны первого вертикала, ПЗ-90.11 П1.2
    XYZ_m[0] = (N+BLH[2])*std::cos(BLH[0])*std::cos(BLH[1]);
    XYZ_m[1] = (N+BLH[2])*std::cos(BLH[0])*std::sin(BLH[1]);
    XYZ_m[2] = ((1-e2)*N + BLH[2])*std::sin(BLH[0]);
    return;
}


/**
 * @brief Функция для преобразования геодезических координат широты (B),
 * долготы (L) и высоты (H) в геоцентрические координаты
 * Используются формулы П1.1 из ПЗ-90.11
 */
void XYZtoBLH(
    const double X_m,
    const double Y_m,
    const double Z_m,
    double & B_rad,
    double & L_rad,
    double & H_m)
{
    return;
}


/**
 * @brief Построение пары кривых зависимости плотности от высоты по параметрам,
 * указанным поользователем в графическом интерфейсе.
 * Геоцентрические координаты X, Y и Z рассчитываются на основе широты и долготы,
 * указанной пользователем, по формулам
 */
void MainWindow::PlotByParameters()
{
    m_chart1->removeAllSeries();

    // QList<QGroupBox*> curveBoxes = {ui->gbCurve1, ui->gbCurve2};

    if(ui->gbCurve1->isChecked() && ui->chkGOSTon->isChecked()){

        QLineSeries * series = new QLineSeries();
        series->setName("Кривая 1");
        double h_km = 120;
        double h_step = 20;
        while(h_km <= 1500){
            double BLH[3] = {
                qDegreesToRadians(ui->spnB_1->value()),
                qDegreesToRadians(ui->spnL_1->value()),
                1000.0 * h_km
            };
            double XYZ[3] = {0.0};
            BLHtoXYZ(BLH, XYZ);
            XYZ[0] = XYZ[0]/1000.0;
            XYZ[1] = XYZ[1]/1000.0;
            XYZ[2] = XYZ[2]/1000.0;

            double rho = atmosGOST_R_25645_166_2004(
                h_km,
                ui->spnF107_1->value(),
                ui->spnKp_1->value(),
                ui->spnF81_1->value(),
                ui->spnDoY_1->value(),
                XYZ,
                ui->spnTs_1->value(),
                ui->spnS_1->value(),
                ui->spnAlpha_1->value(),
                ui->spnDelta_1->value());
            series->append(h_km, rho);
            h_km += h_step;
        }
        m_chart1->addSeries(series);
        series->attachAxis(m_axisY);
        series->attachAxis(m_axisX);
    }

    if(ui->gbCurve2->isChecked() && ui->chkGOSTon->isChecked()){

        QLineSeries * series = new QLineSeries();
        series->setName("Кривая 2");
        double h_km = 120;
        double h_step = 20;
        while(h_km <= 1500){
            double BLH[3] = {
                qDegreesToRadians(ui->spnB_2->value()),
                qDegreesToRadians(ui->spnL_2->value()),
                1000.0 * h_km
            };
            double XYZ[3] = {0.0};
            BLHtoXYZ(BLH, XYZ);
            XYZ[0] = XYZ[0]/1000.0;
            XYZ[1] = XYZ[1]/1000.0;
            XYZ[2] = XYZ[2]/1000.0;
            double rho = atmosGOST_R_25645_166_2004(
                h_km,
                ui->spnF107_2->value(),
                ui->spnKp_2->value(),
                ui->spnF81_2->value(),
                ui->spnDoY_2->value(),
                XYZ,
                ui->spnTs_2->value(),
                ui->spnS_2->value(),
                ui->spnAlpha_2->value(),
                ui->spnDelta_2->value());
            series->append(h_km, rho);
            h_km += h_step;
        }
        m_chart1->addSeries(series);
        series->attachAxis(m_axisY);
        series->attachAxis(m_axisX);
    }

    if(ui->gbCurve1->isChecked() && ui->chkNrlmsise00on->isChecked()){
        QLineSeries * series_nrlmsise00_1 = new QLineSeries();
        series_nrlmsise00_1->setName("Кривая 1 (NRLMSISE-00)");

        double h_km = 120;
        double h_step = 20;

        nrlmsise_input input;
        input.doy = static_cast<int>(ui->spnDoY_1->value());
        input.year=0; /* without effect */
        input.sec=ui->spnTs_1->value();
        input.g_lat=ui->spnB_1->value();
        input.g_long=ui->spnL_1->value();

        // см. комментарий к исходнику структуры nrlmsise_input:
        input.lst=input.sec/3600.0 + input.g_long/15.0;

        input.f107A=ui->spnF81_1->value();
        input.f107=ui->spnF107_1->value();
        input.ap=KpToApLinearInterp( ui->spnKp_1->value());
        qDebug() << "ap = " << input.ap;

        ap_array ap_array;
        for(int i = 0; i < 7; i++)
            ap_array.a[i] = KpToApLinearInterp( ui->spnKp_1->value());
        input.ap_a = &ap_array;

        // по заполнению см. комментарий к исходнику nrlmsise_flags
        nrlmsise_flags flags;
        //flags.switches[0] = 0.0;
        for(int i = 0; i < 24; i++)
            flags.switches[i] = 1.0;
        flags.switches[9] = -1;

        while(h_km <= 1000){
            input.alt=h_km;
            nrlmsise_output output;
            gtd7(&input, &flags, &output);

            double rho = output.d[5];
            series_nrlmsise00_1->append(h_km, rho);
            h_km += h_step;
        }

        m_chart1->addSeries(series_nrlmsise00_1);
        series_nrlmsise00_1->attachAxis(m_axisY);
        series_nrlmsise00_1->attachAxis(m_axisX);
    }

    if(ui->gbCurve2->isChecked() && ui->chkNrlmsise00on->isChecked()){
        QLineSeries * series_nrlmsise00_2 = new QLineSeries();
        series_nrlmsise00_2->setName("Кривая 2 (NRLMSISE-00)");

        double h_km = 120;
        double h_step = 20;

        nrlmsise_input input;
        input.doy = static_cast<int>(ui->spnDoY_2->value());
        input.year=0; /* without effect */
        input.sec=ui->spnTs_2->value();
        input.g_lat=ui->spnB_2->value();
        input.g_long=ui->spnL_2->value();

        // см. комментарий к исходнику структуры nrlmsise_input:
        input.lst=input.sec/3600.0 + input.g_long/15.0;

        input.f107A=ui->spnF81_2->value();
        input.f107=ui->spnF107_2->value();
        input.ap=KpToApLinearInterp( ui->spnKp_2->value());
        qDebug() << "ap = " << input.ap;

        ap_array ap_array;
        for(int i = 0; i < 7; i++)
            ap_array.a[i] = KpToApLinearInterp( ui->spnKp_2->value());
        input.ap_a = &ap_array;

        // по заполнению см. комментарий к исходнику nrlmsise_flags
        nrlmsise_flags flags;
        //flags.switches[0] = 0.0;
        for(int i = 0; i < 24; i++)
            flags.switches[i] = 1.0;
        flags.switches[9] = -1;

        while(h_km <= 1000){
            input.alt=h_km;
            nrlmsise_output output;
            gtd7(&input, &flags, &output);

            double rho = output.d[5];
            series_nrlmsise00_2->append(h_km, rho);
            h_km += h_step;
        }

        m_chart1->addSeries(series_nrlmsise00_2);
        series_nrlmsise00_2->attachAxis(m_axisY);
        series_nrlmsise00_2->attachAxis(m_axisX);
    }
}

