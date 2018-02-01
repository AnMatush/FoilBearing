using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace Charp_mesh_apdl
{
    class Program
    {
        static void Main(string[] args)
        {
            Console.WriteLine("Начало работы");
            //время выполнения программы/////
            DateTime time11 = DateTime.Now;
            for (int i = 0; i < 20000000; i++) { }
            ////////////////////////////////////
            ///////////////////////////////////////////////////////////////////////////////////////////////
            /////////////////////Эти данные потом считать из файла////////////////////////////////////////
            double Shaft_radius = 32;
            double Shaft_radius_excentrisitet = 0;
            double Theta = 0; //угол эксцентриситета 
            //double Delt_Excentrisitet = 0.004;
            double popravka = 0.1;
            //int max_iter = 1;             //максимальное число итераций
            //double deltmax = 0.0005;       //максимальная невязка по зазору
            //double razmer = 1;              //размер конечного элемента
            double dlinna = 60;              //длина подшипника
            int Ndel = 1;                  //разряжённость сетки давлений
            double Sravnenie_y = 0.2;      //сравнивается соотношение соседних узлов yб чтоб присвоить им 1 номер
            double Perepad = 0.1 - 0.01;   //Перепад высот в нахлёсте лепестка
            int num_foil = 4;                //Число лепестков
            double Hmid = 0.00001;      //Характеристика зазора (отчего именно такая?)

            //////////////////////////////////////////////////////////////////////////////////////////////



            /*Чтение из файла*/
            FileStream file1 = new FileStream("D:\\bearing\\BBEEAARRIINNGG\\Number_node_coord1.txt", FileMode.Open, FileAccess.Read); //открывает файл только на чтение
            StreamReader reader = new StreamReader(file1, System.Text.Encoding.Default); // создаем «потоковый читатель» и связываем его с файловым потоком 
            string[] AllDataS = reader.ReadToEnd().Split(new char[] { ' ' }, StringSplitOptions.RemoveEmptyEntries); //создание одномерного массива символов 

            int MaxNMass = AllDataS.Length; //максимальное число элементов массива
            int Num_strok = AllDataS.Length / 5; //Число строе в массиве


            //
            //Создание нескольких массивов по столбцам
            //
            int a = 0;
            //Инициация массивов
            double[] Node_Number = new double[Num_strok];
            double[] Coord_x = new double[Num_strok];
            double[] Coord_y = new double[Num_strok];
            double[] Coord_z = new double[Num_strok];
            double[] Deform = new double[Num_strok];
            double[] New_x = new double[Num_strok];
            double[] Zazor = new double[Num_strok];
            double[,] All_coord_number = new double[Num_strok, 7];

            //Создание массивов по столбцам

            //Массив номеров узлов
            for (int i = 0; i < MaxNMass; i = i + 5)
            {
                Node_Number[a] = Math.Round(Convert.ToDouble(AllDataS[i], System.Globalization.CultureInfo.InvariantCulture));
                a++;
            }
            a = 0;
            //Массив координат x
            for (int i = 1; i < MaxNMass; i = i + 5)
            {
                Coord_x[a] = Convert.ToDouble(AllDataS[i], System.Globalization.CultureInfo.InvariantCulture);
                a++;
            }
            a = 0;

            //Массив координат y c округлением до целого
            for (int i = 2; i < MaxNMass; i = i + 5)
            {
                Coord_y[a] = Convert.ToDouble(AllDataS[i], System.Globalization.CultureInfo.InvariantCulture);
                a++;
            }
            a = 0;

            //Массив координат z
            for (int i = 3; i < MaxNMass; i = i + 5)
            {
                Coord_z[a] = Convert.ToDouble(AllDataS[i], System.Globalization.CultureInfo.InvariantCulture);
                a++;
            }
            a = 0;

            //Массив деформаций
            for (int i = 4; i < MaxNMass; i = i + 5)
            {
                Deform[a] = Convert.ToDouble(AllDataS[i], System.Globalization.CultureInfo.InvariantCulture);
                a++;
            }

            //Массив координат х вала с эксцентриситетом
            for (int i = 0; i < Num_strok; i++)
            {
                New_x[i] = Shaft_radius + Shaft_radius_excentrisitet * Math.Cos((Coord_y[i] + Theta) * Math.PI / 180);
            }

            //Массив Зазоров
            for (int i = 0; i < Num_strok; i++)
            {
                //Zazor[i] = Coord_x[i] + Deform[i] - New_x[i] - popravka;
                Zazor[i] = (Coord_x[i] - New_x[i] - popravka);
            }

            //Создание итогового массива №, x, y, z, x+excentr, +
            for (int i = 0; i < Num_strok; i++)
            {
                for (int j = 0; j <= 6; j++)
                {
                    if (j == 0) All_coord_number[i, j] = Node_Number[i];
                    if (j == 1) All_coord_number[i, j] = Coord_x[i];
                    if (j == 2) All_coord_number[i, j] = Coord_y[i];
                    if (j == 3) All_coord_number[i, j] = Coord_z[i];
                    if (j == 4) All_coord_number[i, j] = Deform[i];
                    if (j == 5) All_coord_number[i, j] = New_x[i]; //Эксцентриситет+x
                    if (j == 6) All_coord_number[i, j] = Zazor[i]/1000; //Зазоры
                }
            }


            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////          
            //Сортировка массива All_coord_number методом расчёски по увеличению координаты y
            double fakt = 1.2473309; //Фактор уменьшения
            double gap = Num_strok; //Длинна массива
            bool obmen = true;   //Флаг
            while (gap > 1 || obmen)
            {
                gap = gap / fakt;
                if (gap < 1) { gap = 1; }
                int i = 0;
                obmen = false;
                while (i + gap < Num_strok)
                {
                    int igap = i + (int)gap;
                    if (Math.Abs(All_coord_number[i, 2] - All_coord_number[igap, 2]) < Sravnenie_y) All_coord_number[i, 2] = All_coord_number[igap, 2];
                    if (All_coord_number[i, 2] > All_coord_number[igap, 2])
                    {
                        for (int j = 0; j <= 6; j++)
                        {
                            double[] oben_in = new double[7];
                            oben_in[j] = All_coord_number[i, j];
                            All_coord_number[i, j] = All_coord_number[igap, j];
                            All_coord_number[igap, j] = oben_in[j];
                            obmen = true;
                        }
                    }
                    //Сортировка оставшихся элементов по увеличению Z
                    if (All_coord_number[i, 2] == All_coord_number[igap, 2] & All_coord_number[i, 3] > All_coord_number[igap, 3])
                    {
                        for (int j = 0; j <= 6; j++)
                        {
                            double[] oben_in = new double[7];
                            oben_in[j] = All_coord_number[i, j];
                            All_coord_number[i, j] = All_coord_number[igap, j];
                            All_coord_number[igap, j] = oben_in[j];
                            obmen = true;
                        }
                    }
                    i++;
                }
            }






            /////////////////ИСПРАВИТЬ ЧТО ТО НЕ РАБОТАЕТ////////////////
            /////Подсчёт числа элементов по Y
            int max_iter_y = 0;
            a = 0;
            for (int i = 1; i < Num_strok; i++)
            {
                if (Math.Abs(All_coord_number[i, 2] - All_coord_number[i - 1, 2]) > Sravnenie_y) max_iter_y = a = a + 1;
            }
            max_iter_y = max_iter_y - 5;
            double delta_y = 360 / max_iter_y; //шаг между соседними узлами в одной окружности


            /////Подсчёт числа элементов по Z
            a = 1;
            int max_iter_z = 0;
            while (a < max_iter_y & All_coord_number[a, 2] == All_coord_number[a - 1, 2]) { max_iter_z = a = a + 1; }

            ///////Подсчёт числа стыков лепестков и номера массива All_coord_number с которых начинаются новые элементы
            int max_iter_x = 0;
            a = 0;
            int[] num_i_of_zazor = new int[num_foil];//Одномерный смассив номеров All_coord_number с которых начинается новый лепесток
            for (int i = 1; i < Num_strok; i++)
                if (Math.Abs(All_coord_number[i, 1] - All_coord_number[i - 1, 1]) > Perepad) { max_iter_x = max_iter_x + 1; num_i_of_zazor[a] = i + 1; a++; }
            int max_num_foil = num_i_of_zazor[2] - num_i_of_zazor[1];
            int num_y_of_foil = max_iter_y / num_foil;


            //Преобразование из массива All_coord_number формата kх6
            //В массив H[i, j], где i-массив углов y, j-массив длинны z
            //Осуществляется по формуле k=max_iter_z*i+j
            //Где max_iter_z - число элементов в массиве j для 1 угла y 
            // for (int i = 0; i < ny; i++)
            // {
            //    for (int j = 0; j < nz; j++)
            //    {
            //        PXZ[i, j] = 1.0;
            //        H[i, j] = All_coord_number[max_iter_z * i + j, 6];
            //    }
            // }


            //Преобразование из массива All_coord_number формата kх6
            //В массив зазоров H[i, j], где i-массив углов y, j-массив длинны z
            //Осуществляется по формуле k=max_iter_z*i+j
            //Где max_iter_z - число элементов в массиве j для 1 угла y 
            //Создание массивов для первых лепестков
            double[,,] H = new double[num_y_of_foil, max_iter_z, num_foil];     //Массив зазоров
            double[,,] PXZ = new double[num_y_of_foil, max_iter_z, num_foil];   //Массив невязок


            for (int k = 1; k < num_foil; k++)
            {
                for (int i = 0; i < num_y_of_foil; i++)
                {
                    for (int j = 0; j < max_iter_z; j++)
                    {
                        H[i, j, k - 1] = All_coord_number[(max_iter_z * i + j) + num_i_of_zazor[k - 1], 6];
                        PXZ[i, j, k - 1] = 1;
                    }
                }
            }
            //Создание массива для последнего лепестка (переходящего через конец массива All_coord_number числа узлов )
            for (int i = 0; i < num_y_of_foil; i++)
            {
                for (int j = 0; j < max_iter_z; j++)
                {
                    if ((max_iter_z * i + j) + num_i_of_zazor[num_foil - 1] < Num_strok) { H[i, j, num_foil - 1] = All_coord_number[(max_iter_z * i + j) + num_i_of_zazor[num_foil - 1], 6]; PXZ[i, j, num_foil - 1] = 1; }
                    if ((max_iter_z * i + j) + num_i_of_zazor[num_foil - 1] > Num_strok)
                    {
                        int istart = 0;
                        H[i, j, num_foil - 1] = All_coord_number[(max_iter_z * istart + j), 6]; PXZ[i, j, num_foil - 1] = 1;
                        istart++;
                    }
                }
            }

            for (int i=0; i<50; i++)
            {
                for (int j = 0; j < 61; j++)
                    H[i, j, 0] = (0.01/1000 - (0.01/1000 - 0.001/1000) / 50 * i)/0.000002;
            }

            for (int i = 0; i < num_y_of_foil; i++) { Console.WriteLine("{0} , {1} , {2} , {3} ", H[i, 20, 0], H[i, 20, 1], H[i, 20, 2], H[i, 20, 3]); }



            Console.WriteLine("Начало расчёта уравнений рейнольдса");
            //Сбор и рассчёт массива давлений
            double[,,] Pressure = new double[num_y_of_foil, max_iter_z,num_foil];
                for (int k = 0; k < num_foil; k++)
                    {
                        Console.WriteLine("Решение для лепестка: " + (k+1));
                        Pressure = Pressure_solve(num_y_of_foil, max_iter_z, Shaft_radius, dlinna, num_foil, Ndel, delta_y, H, PXZ, k);
                    }
            for (int i = 0; i < num_y_of_foil; i++) { Console.WriteLine("{0} , {1} , {2} , {3} ", Pressure[i, 20, 0], Pressure[i, 20, 1], Pressure[i, 20, 2], Pressure[i, 20, 3]); }

            double SUM = 0;
            for (int i = 0; i < 50; i++)
            {
                SUM = Pressure[i, 10, 0] + SUM;
            }
            Console.WriteLine("Сумма давлений " + SUM);


            //Преобразование из массива Pressure формата ny х nz
            //В массив Massiv_Pressure[i, 2], где i-число узлов, 1-номера узлов, 2-давления
            //Осуществляется по формуле для ny массива PPXZ деление с округлением до целого вниз i/max_iter_z
            // для nz массива PPXZ i-(деление с округлением до целого вниз i/ max_iter_z)*max_iter_z
            //Где max_iter_z - число элементов в массиве j для 1 угла y 

            //сбор для лепестков кроме последнего
            double[,] Massiv_Pressure = new double[Num_strok, 2];
             for (int k = 0; k < num_foil-1 ; k++)
             {
               int j = 0;
               for (int i = num_i_of_zazor[k]; i < num_i_of_zazor[k+1]; i++)
                 {
                    if (j== max_num_foil) { j = 0; }
                 Massiv_Pressure[i, 0] = All_coord_number[i, 0];
                 Massiv_Pressure[i, 1] = Pressure[j / max_iter_z, j - (j / max_iter_z) * max_iter_z,k];
                    j++;
                }
                
             }
            //сбор для последнего лепестка
            int num = num_foil - 1;
                int count = 0;
                for (int i = num_i_of_zazor[num]; i < Num_strok; i++)
                {
                    Massiv_Pressure[i, 0] = All_coord_number[i, 0];
                    Massiv_Pressure[i, 1] = Pressure[count / max_iter_z, count - (count / max_iter_z) * max_iter_z, num];
                    count++;
                }
                
                for (int i = 0; i < num_i_of_zazor[num]; i++)
                {
                    if (count == max_num_foil) { count = 0; }
                    Massiv_Pressure[i, 0] = All_coord_number[i, 0];
                    Massiv_Pressure[i, 1] = Pressure[count / max_iter_z, count - (count / max_iter_z) * max_iter_z, num];
                count++;
                }
           


           // for (int i = num_i_of_zazor[2]; i < num_i_of_zazor[3]; i++) { Console.WriteLine(Massiv_Pressure[i, 1]); }
          
            //////////////////////////////////////////////////////////////
            ////////////////////Создание массива давлений для ANSYS//////
            ///////////////////И его вывод в файл////////////////////////
            //////////////////////////////////////////////////////////////
            //////////////////////////////////////////////////////////////

          
                        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////          
                        //Сортировка массива Massiv_Pressure методом расчёски по увеличению номеров узлов
                        double fakt_1 = 1.2473309; //Фактор уменьшения
                        double gap_1 = Num_strok; //Длинна массива
                        bool obmen_1 = true;   //Флаг
                        while (gap_1 > 1 || obmen_1)
                        {
                            gap_1 = gap_1 / fakt_1;
                            if (gap_1 < 1) { gap_1 = 1; }
                            int i = 0;
                            obmen_1 = false;
                            while (i + gap_1 < Num_strok)
                            {
                                int igap_1 = i + (int)gap_1;
                                if (Massiv_Pressure[i, 0] > Massiv_Pressure[igap_1, 0])
                                {
                                    for (int j = 0; j <= 1; j++)
                                    {
                                        double[] oben_in = new double[7];
                                        oben_in[j] = Massiv_Pressure[i, j];
                                        Massiv_Pressure[i, j] = Massiv_Pressure[igap_1, j];
                                        Massiv_Pressure[igap_1, j] = oben_in[j];
                                        obmen_1 = true;
                                    }
                                }
                                i++;
                            }
                        }
          


                        FileStream file_pressure = new FileStream("D:\\bearing\\BBEEAARRIINNGG\\Pressure_CSharp.txt", FileMode.Create, FileAccess.Write); //открывает файл на запись
                        StreamWriter writer_pressure = new StreamWriter(file_pressure, System.Text.Encoding.Default); // создаем «потоковый читатель» и связываем его с файловым потоком 
                        for (int i = 0; i < Num_strok; i++)
                        {
                            System.Threading.Thread.CurrentThread.CurrentCulture = new System.Globalization.CultureInfo("en-US");
                            writer_pressure.WriteLine((Massiv_Pressure[i, 1]));
                        }
           
            /*
            FileStream test = new FileStream("D:\\bearing\\BBEEAARRIINNGG\\TTesttt.txt", FileMode.Create, FileAccess.Write); //открывает файл на запись
            StreamWriter writer_test = new StreamWriter(test, System.Text.Encoding.Default); // создаем «потоковый читатель» и связываем его с файловым потоком 
            for (int i = 0; i < ny; i++)
            {
                writer_test.WriteLine("{0},{1}", H[i, 30], PPXZ[i, 30]);
            }
            */


            Console.WriteLine(max_iter_y);
            //Console.WriteLine(max_iter_z);
            //Console.WriteLine(max_iter_x);
            //for (int i=0; i< num_y_of_foil; i++) { Console.WriteLine(Pressure[i,20]); }

            DateTime time22 = DateTime.Now;
            Console.WriteLine("Время выполнения суммарное: {0}", (time22 - time11).Seconds);

            Console.ReadLine();
            /*
                        

                     
                     
                        


                       
                      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
                      ////////////////////////// Вывод в файл массива отсортированного
                      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
                      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
                      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 

                       FileStream file_out = new FileStream("D:\\bearing\\BBEEAARRIINNGG\\Number_node_coord_CSharp.txt", FileMode.Create, FileAccess.Write); //открывает файл только на чтение
                      StreamWriter writer = new StreamWriter(file_out, System.Text.Encoding.Default); // создаем «потоковый читатель» и связываем его с файловым потоком 
                                  for (int i = 0; i < Num_strok; i++)
                                  {

                                          writer.WriteLine("{0} | {1} | {2} | {3} | {4} | {5} | {6}", All_coord_number[i, 0], All_coord_number[i, 1], All_coord_number[i, 2], All_coord_number[i, 3],All_coord_number[i, 4], All_coord_number[i, 5], All_coord_number[i, 6]);

                                  }
                       ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
                       ////////////////////////// Вывод в файл массива отсортированного
                      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
                      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
                      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    */
        }
        ///////////////////Подпрограмма 1//////////////////////
        public static double [,,] Pressure_solve (int max_iter_y, int max_iter_z, double Shaft_radius, double dlinna, int num_of_foil, int Ndel, double delta_y, double[,,] Hg, double [,,]PXZg, int k)
        {
            int ny = max_iter_y;                //число узлов по y - 4(по числу пересечений лепестков)?????????????????????????????????????
            int nz = max_iter_z;                //число узлов по z
            double Shaft = Shaft_radius / 1000;        //Радиус вала (м)
            double BLep = dlinna / 1000;      //Длинна подшипника
            int z = num_of_foil;             //число лепестков

            double alfa = 0.8;            //Коэффициент релаксации
            int Nit = 3;                //Число итераций прогона
            int Nit1 = 800;                     //число итераций расчёта давлений
            int Pa = 100000;            //Атмосферное давление 1атм
            int Omega = 3000;           //Частота вращения 
           
            double Myu = 0.0000183;     //Вязкость
            double Cz = 120;            //Константа Сазерленда
            double TA = 300;            //Температура смазочного слоя в подшипнике (К)
            double Hmid = 0.002/1000;      //Характеристика зазора (отчего именно такая?)

            double Myu0;                //Вязкость с поправкой на температуру
            double dx;                  //Угловой шаг узлов в подшипнике (радиан)
            double DZ;                  //Относительный шаг элементов по оси Z
            double dtau;                //время релаксации 
            double lamb;                //Коэффициент сжимаемости
            double OmegaK = 100;        //частота осциляций 3000
            double SIGMA;               //Коэффициент сдавливания

            Myu0 = Myu * (273.0 + Cz) / ((TA) + Cz) * Math.Pow((TA) / 273.0, 1.5);
            Console.WriteLine("Мю0 " + Myu0);
            dx = delta_y*2*Math.PI/180;
            Console.WriteLine("Угловой шаг узлов в подшипнике (радиан) " + dx);
            DZ = BLep / Shaft / (nz - 1) * Ndel;
            Console.WriteLine("Относительный шаг элементов по оси Z " + DZ);
            dtau = 1 / (2 * 1 * (1 / Math.Pow(dx, 2) + 1 / Math.Pow(DZ, 2)));
            Console.WriteLine("время релаксации  " + dtau);
            lamb = 6 * Myu * Omega * Math.Pow(Shaft, 2) / Pa / Hmid / Hmid;
            Console.WriteLine("Число сжимаемости " + lamb);
            SIGMA = 12 * Myu * Math.Pow(Shaft, 2) * OmegaK / Pa / Hmid / Hmid;
            Console.WriteLine("Коэффициент сдавливания " + SIGMA);

            double[,,] PPXZ = new double[ny, nz,z];
            
            for (int it = 1; it <= Nit; it++)
            {
                PPXZ = ReynoldsUstanovlenie(dx, DZ, lamb, dtau, SIGMA, alfa, Nit1, ny, nz, Hg, PXZg, k, num_of_foil);
                Console.WriteLine("прогон: " + it);
            }
            return PPXZ;
        }



        //////////////////Подпрограмма 2//////////////////////////
        public static double[,,] ReynoldsUstanovlenie(double dx, double DZ, double lamb, double dtau, double SIGMA, double alfa, double Nit1, int ny, int nz, double[,,] H, double[,,] PXZ, int k,int num_of_foil)
        {
           //Переменные
            double u1;
            double u2;
            double u3;
            double u4;
            double u5;
            double u6;
            double Uxp1;
            double Uxm1;
            double Uzm1;
            double Uzp1;
            double DelU;
            double PN;
            
            //Объявление зазоров
            double[,,] H3p = new double[ny, nz,num_of_foil];
            double[,,] H3m = new double[ny, nz,num_of_foil];

            for (int it1 = 0; it1 <= Nit1; it1++)
            {
                for (int i = 1; i < ny-1; i++)
                {
                   
                    for (int j = 1; j < nz - 1; j++)
                    {
                        H3p[i, j,k] = Math.Pow((((H[i, j,k] + H[i+1, j,k])) / 2), 3);
                        H3m[i, j,k] = Math.Pow((((H[i, j,k] + H[i-1, j,k])) / 2), 3);
                    }
                }

                for (int i1 = 1; i1 < ny-1 ; i1++) //условия открытого по оси подшипника
                {
                    int I = i1;
                    int im1 = I - 1;
                    int ip1 = I + 1;
                   
                    // if (i1==1) { im1 = 1; }
                    // if(i1==ny-1) { ip1 = ny; }

                    for (int j = 1; j < nz - 1; j++)
                    {
                          
                        Uxp1 = H3p[I, j,k] * (PXZ[ip1, j,k] + PXZ[I, j,k]) / 2 * (PXZ[ip1, j,k] - PXZ[I, j,k]) / dx;
                        u1 = Uxp1;
                        Uxp1 = Uxp1 - lamb * (PXZ[ip1, j,k] + PXZ[I, j,k]) / 2 * (H[ip1, j,k] + H[I, j,k]) / 2;
                        u2 = Uxp1;
                        Uxm1 = H3m[I, j,k] * (PXZ[im1, j,k] + PXZ[I, j,k]) / 2 * (-PXZ[im1, j,k] + PXZ[I, j,k]) / dx;
                        u3 = Uxm1;
                        Uxm1 = Uxm1 - lamb * (PXZ[im1, j,k] + PXZ[I, j,k]) / 2 * (H[im1, j,k] + H[I, j,k]) / 2;
                        u4 = Uxm1;
                        Uzp1 = Math.Pow(((H[I, j + 1,k] + H[I, j,k]) / 2), 3) * (PXZ[I, j + 1,k] + PXZ[I, j,k]) / 2 * (PXZ[I, j + 1,k] - PXZ[I, j,k]) / (DZ);
                        u5 = Uzp1;
                        Uzm1 = Math.Pow(((H[I, j - 1,k] + H[I, j,k]) / 2), 3) * (PXZ[I, j - 1,k] + PXZ[I, j,k]) / 2 * (-PXZ[I, j - 1,k] + PXZ[I, j,k]) / (DZ);
                        u6 = Uzp1;
                        DelU = (Uxp1 - Uxm1) / dx + (Uzp1 - Uzm1) / (DZ);

                        PN = 1 / H[I, j,k] * (DelU * dtau / SIGMA + H[I, j,k] * PXZ[I, j,k]);
                        PXZ[I, j,k] = alfa * PN + (1 - alfa) * PXZ[I, j,k];
                    }
                }
            }
            return PXZ;
         }

    }
}
