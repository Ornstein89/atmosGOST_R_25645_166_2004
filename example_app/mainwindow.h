#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QChart>
#include <QLogValueAxis>
#include <QMainWindow>
#include <QValueAxis>

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

public slots:
    void PlotBunchOfCurves();
    void PlotByParameters();

    void SetupChart();

    void SetupGUI();

private:
    Ui::MainWindow *ui;

    QChart * m_chart1 = nullptr;

    QValueAxis * m_axisX = nullptr;
    QLogValueAxis * m_axisY = nullptr;

};
#endif // MAINWINDOW_H
