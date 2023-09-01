#include "mainwindow.h"
#include "./ui_mainwindow.h"

#include <QChart>
#include <QChartView>
#include <QLineSeries>
#include <QLogValueAxis>
#include <QValueAxis>

#include <atmosGOST_R_25645_166_2004.h>



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

void MainWindow::PlotByParameters()
{
    m_chart1->removeAllSeries();

    // QList<QGroupBox*> curveBoxes = {ui->gbCurve1, ui->gbCurve2};

    if(ui->gbCurve1->isChecked()){

        QLineSeries * series = new QLineSeries();
        series->setName("Кривая 1");
        double h_km = 120;
        double h_step = 20;
        while(h_km <= 1500){

            double X[3] = {
                           ui->spnX_1->value(),
                           ui->spnY_1->value(),
                           ui->spnZ_1->value()};
            double rho = atmosGOST_R_25645_166_2004(
                h_km,
                ui->spnF107_1->value(),
                ui->spnKp_1->value(),
                ui->spnF81_1->value(),
                ui->spnDoY_1->value(),
                X,
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

    if(ui->gbCurve1->isChecked()){

        QLineSeries * series = new QLineSeries();
        series->setName("Кривая 2");
        double h_km = 120;
        double h_step = 20;
        while(h_km <= 1500){
            double X[3] = {
                           ui->spnX_2->value(),
                           ui->spnY_2->value(),
                           ui->spnZ_2->value()};
            double rho = atmosGOST_R_25645_166_2004(
                h_km,
                ui->spnF107_2->value(),
                ui->spnKp_2->value(),
                ui->spnF81_2->value(),
                ui->spnDoY_2->value(),
                X,
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
}

