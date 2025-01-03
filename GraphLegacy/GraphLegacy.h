#pragma once
#include <iostream>
#include <vector>
#include <string>
#include <limits>
#include <queue>
#include <functional>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cassert>
#include <random>
/*
ИСПРАВЛЕНО:
1. КОММЕНТАРИИ: ДОБАВЛЕНО --- ОПИСАНИЕ БРОСАНИЯ ИСКЛЮЧЕНИЙ
2. КОД: РАНЬШЕ ГРАФ ХРАНИЛ ПОЛЕ РАЗМЕРА. СЕЙЧАС -- НЕ ХРАНИТ. ДРОБАВИЛ CONST, ГДЕ НЕОБХОДИМО
3. МИН.ПУТЬ: ЛИСТОЧЕК С ОПИСАНИЕМ АЛГОРИТМА

АНАЛИЗ ГРАФА:
1. ОПИСАНИЕ --- ДОБАВЛЕНО ОПИСАНИЕ ТОГО, ЧТО ЗНАЧАТ СВЯЗИ МЕЖДУ УЗЛАМИ (УЧЕНЫМИ)
2. МЕТРИКИ --- ОПИСАНИЕ ТОГО, ЧТО ТАКОЕ РАЗНЫЕ ПАРАМЕТРЫ В ЛАБОРАТОРИИ ДАННЫХ.

*/



/*
| Операция                                   | Сложность BigO         |
|--------------------------------------------|-------------------------|
| findVertex(const T& vertex)                | O(n)                    |
| getVertexPosition(const T& vertex)         | O(n)                    |
| insertVertex(const T& vertex)              | O(n^2)                  |
| removeVertex(const T& vertex)              | O(n^2)                  |
| insertEdge(const T& vertex1, const T& vertex2, double weight) | O(n) |
| removeEdge(const T& vertex1, const T& vertex2) | O(n)                |
| getAdjacentVertices(const T& vertex)       | O(n)                    |
| getEdgeWeight(const T& vertex1, const T& vertex2) | O(n)               |
| printVertices()                            | O(n)                    |
| printEdges()                               | O(n^2)                  |
| getEdgeCount()                             | O(n^2)                  |
| getVertexCount()                           | O(1)                    |
| findShortestPath(const T& startVertex, const T& endVertex) | O(log(n + E)) |
| depthFirstSearch(const T& startVertex)     | O(n + E)                |
| breadthFirstSearch(const T& startVertex)   | O(n + E)                |
| printAdjacencyMatrix()                     | O(n^2)                  |
| updateEdge(const T& startVertex, const T& endVertex, double weight) | O(n) |
| findAllPaths(const T& startVertex, const T& endVertex) | O(n! + E)        |
| getVertexDegree(const T& vertex)           | O(n)                    |
| getEdgesWeightSum()                        | O(n^2)                  |
*/
using namespace std;
// Класс Взвешенного Графа
template <typename T>
class WeightedGraphLegacy {
private:
    vector<T> vertices; // Список вершин
    vector<vector<double>> adjacencyMatrix; // Матрица смежности. Хранится нижняя треугольная матрица для экономии памяти

public:
    // Конструктор
    WeightedGraphLegacy() {}

    // Деструктор
    ~WeightedGraphLegacy() {}

    // Метод поиска вершины в списке
    int findVertex(const T& vertex) const {
        for (int i = 0; i < vertices.size(); ++i) {
            if (vertices[i] == vertex) {
                return i;
            }
        }
        return -1; // Вершина не найдена
    }

    // Метод поиска позиции вершины в списке
    int getVertexPosition(const T& vertex) const {
        return findVertex(vertex);
    }

    // Метод вставки вершины. Бросает исключение out_of_range если вершина существует 
    void insertVertex(const T& vertex) {
        if (findVertex(vertex) != -1) {
            throw out_of_range("Вершина уже существует в графе.");
        }

        int newSize = vertices.size() + 1;
        vertices.push_back(vertex);

        // Обновление матрицы смежности
        if (adjacencyMatrix.empty()) {
            adjacencyMatrix.resize(newSize, vector<double>(newSize, 0.0));
        }
        else {
            for (auto& row : adjacencyMatrix) {
                row.resize(newSize, 0.0);
            }
            adjacencyMatrix.resize(newSize, vector<double>(newSize, 0.0));
        }
    }

    // Метод удаления вершины. Бросает исключение out_of_range если вершина не существует 
    void removeVertex(const T& vertex) {
        int position = findVertex(vertex);
        if (position == -1) {
            throw out_of_range("Вершина не найдена в графе.");
        }

        vertices.erase(vertices.begin() + position);

        // Обновление матрицы смежности
        adjacencyMatrix.erase(adjacencyMatrix.begin() + position);
        for (int i = 0; i < vertices.size(); ++i) {
            adjacencyMatrix[i].erase(adjacencyMatrix[i].begin() + position);
        }
    }

    //// Метод вставки ребра
    //void insertEdge(const T& vertex1, const T& vertex2, double weight) {
    //    int position1 = findVertex(vertex1);
    //    int position2 = findVertex(vertex2);

    //    if (position1 == -1 || position2 == -1) {
    //        throw out_of_range("Одна или обе вершины не найдены в графе.");
    //    }

    //    adjacencyMatrix[position1][position2] = weight;
    //    adjacencyMatrix[position2][position1] = weight; // Для неориентированного графа
    //}

    //// Метод удаления ребра
    //void removeEdge(const T& vertex1, const T& vertex2) {
    //    int position1 = findVertex(vertex1);
    //    int position2 = findVertex(vertex2);

    //    if (position1 == -1 || position2 == -1) {
    //        throw out_of_range("Одна или обе вершины не найдены в графе.");
    //    }

    //    adjacencyMatrix[position1][position2] = 0.0;
    //    adjacencyMatrix[position2][position1] = 0.0; // Для неориентированного графа
    //}

    // Метод вставки ребра. Бросает исключение out_of_range если вершина не существует 
    void insertEdge(const T& vertex1, const T& vertex2, double weight) {
        int position1 = findVertex(vertex1);
        int position2 = findVertex(vertex2);

        if (position1 == -1 || position2 == -1) {
            throw out_of_range("Одна или обе вершины не найдены в графе.");
        }

        if (position1 > position2) {
            swap(position1, position2);
        }

        adjacencyMatrix[position2][position1] = weight;
    }

    // Метод удаления ребра. Бросает исключение out_of_range если вершина не существует 
    void removeEdge(const T& vertex1, const T& vertex2) {
        int position1 = findVertex(vertex1);
        int position2 = findVertex(vertex2);

        if (position1 == -1 || position2 == -1) {
            throw out_of_range("Одна или обе вершины не найдены в графе.");
        }

        if (position1 > position2) {
            swap(position1, position2);
        }

        adjacencyMatrix[position2][position1] = 0.0;
    }


    // Метод возврата списка смежных вершин. Бросает исключение out_of_range если вершина не существует 
    vector<T> getAdjacentVertices(const T& vertex) const {
        int position = findVertex(vertex);
        if (position == -1) {
            throw out_of_range("Одна или обе вершины не найдены в графе.");
        }

        vector<T> adjacentVertices;
        for (int i = 0; i < vertices.size(); ++i) {
            if (adjacencyMatrix[position][i] != 0.0) {
                adjacentVertices.push_back(vertices[i]);
            }
        }
        return adjacentVertices;
    }

    //// Метод возврата веса ребра между двумя вершинами
    //double getEdgeWeight(const T& vertex1, const T& vertex2) {
    //    int position1 = findVertex(vertex1);
    //    int position2 = findVertex(vertex2);

    //    if (position1 == -1 || position2 == -1) {
    //        throw out_of_range("Одна или обе вершины не найдены в графе.");
    //    }

    //    return adjacencyMatrix[position1][position2];
    //}

    // Метод возврата веса ребра между двумя вершинами. Бросает исключение out_of_range если вершина не существует 
    double getEdgeWeight(const T& vertex1, const T& vertex2) const {
        int position1 = findVertex(vertex1);
        int position2 = findVertex(vertex2);

        if (position1 == -1 || position2 == -1) {
            throw out_of_range("Одна или обе вершины не найдены в графе.");
        }

        if (position1 > position2) {
            swap(position1, position2);
        }

        return adjacencyMatrix[position2][position1];
    }

    // Метод вывода списка вершин
    void printVertices() const {
        cout << "Список вершин: ";
        for (const auto& vertex : vertices) {
            cout << vertex << " ";
        }
        cout << endl;
    }

    //// Метод вывода списка ребер
    //void printEdges() {
    //    cout << "Список ребер: ";
    //    for (int i = 0; i < size; ++i) {
    //        for (int j = i + 1; j < size; ++j) {
    //            if (adjacencyMatrix[i][j] != 0.0) {
    //                cout << vertices[i] << " - " << vertices[j] << " (вес: " << adjacencyMatrix[i][j] << ") ";
    //            }
    //        }
    //    }
    //    cout << endl;
    //}

    //// Метод возвращения количества ребер
    //int getEdgeCount() {
    //    int edgeCount = 0;
    //    for (int i = 0; i < size; ++i) {
    //        for (int j = i + 1; j < size; ++j) {
    //            if (adjacencyMatrix[i][j] != 0.0) {
    //                edgeCount++;
    //            }
    //        }
    //    }
    //    return edgeCount;
    //}

    // Метод вывода списка ребер
    void printEdges() const {
        cout << "Список ребер: ";
        for (int i = 0; i < vertices.size(); ++i) {
            for (int j = 0; j < i; ++j) {
                if (adjacencyMatrix[i][j] != 0.0) {
                    cout << vertices[i] << " - " << vertices[j] << " (вес: " << adjacencyMatrix[i][j] << ") ";
                }
            }
        }
        cout << endl;
    }

    // Метод возвращения количества ребер
    int getEdgeCount() const {
        int edgeCount = 0;
        for (int i = 0; i < vertices.size(); ++i) {
            for (int j = 0; j < i; ++j) {
                if (adjacencyMatrix[i][j] != 0.0) {
                    edgeCount++;
                }
            }
        }
        return edgeCount;
    }

    // Метод возвращения количества вершин
    int getVertexCount() const {
        return vertices.size();
    }

    //// Процедура поиска минимально возможного пути между двумя вершинами методом Дейкстры
    //void findShortestPath(const T& startVertex, const T& endVertex) {
    //    // Находим позиции стартовой и конечной вершин в графе
    //    int startPosition = findVertex(startVertex);
    //    int endPosition = findVertex(endVertex);

    //    // Если одна или обе вершины не найдены, бросаем исключение
    //    if (startPosition == -1 || endPosition == -1) {
    //        throw out_of_range("Одна или обе вершины не найдены в графе.");
    //    }

    //    // Инициализируем расстояния до всех вершин как бесконечность
    //    vector<double> distances(size, numeric_limits<double>::infinity());
    //    // Исключение: расстояние до стартовой вершины равно 0
    //    distances[startPosition] = 0.0;

    //    // Создаем массив для отслеживания посещенных вершин
    //    vector<bool> visited(size, false);

    //    // Повторяем следующие шаги до тех пор, пока не будут посещены все вершины
    //    for (int i = 0; i < size; ++i) {
    //        // Находим вершину с минимальным расстоянием, которая еще не была посещена
    //        int minDistancePosition = -1;
    //        double minDistance = numeric_limits<double>::infinity();
    //        for (int j = 0; j < size; ++j) {
    //            if (!visited[j] && distances[j] < minDistance) {
    //                minDistancePosition = j;
    //                minDistance = distances[j];
    //            }
    //        }

    //        // Если такая вершина не найдена, алгоритм завершается
    //        if (minDistancePosition == -1) {
    //            break;
    //        }

    //        // Помечаем вершину как посещенную
    //        visited[minDistancePosition] = true;

    //        // Для всех соседних вершин, которые еще не были посещены, рассчитываем новое расстояние
    //        for (int j = 0; j < size; ++j) {
    //            if (!visited[j] && adjacencyMatrix[minDistancePosition][j] != 0.0) {
    //                double newDistance = distances[minDistancePosition] + adjacencyMatrix[minDistancePosition][j];
    //                // Если новое расстояние меньше текущего расстояния, обновляем его
    //                if (newDistance < distances[j]) {
    //                    distances[j] = newDistance;
    //                }
    //            }
    //        }
    //    }

    //    // После завершения алгоритма, расстояние до конечной вершины является минимально возможным путем
    //    if (distances[endPosition] == numeric_limits<double>::infinity()) {
    //        throw out_of_range("Путь между вершинами не найден.");
    //    }
    //    else {
    //        cout << "Минимально возможный путь между вершинами: " << distances[endPosition] << endl;
    //    }
    //}

    // Процедура поиска минимально возможного пути между двумя вершинами методом Дейкстры. Бросает исключение out_of_range если вершина не существует 
    vector<double> findShortestPath(const T& startVertex, const T& endVertex) const {
        // Находим позиции стартовой и конечной вершин в графе
        int startPosition = findVertex(startVertex);
        int endPosition = findVertex(endVertex);

        // Если одна или обе вершины не найдены, бросаем исключение
        if (startPosition == -1 || endPosition == -1) {
            throw out_of_range("Одна или обе вершины не найдены в графе.");
        }

        // Инициализируем расстояния до всех вершин как бесконечность
        vector<double> distances(vertices.size(), numeric_limits<double>::infinity());
        // Исключение: расстояние до стартовой вершины равно 0
        distances[startPosition] = 0.0;

        // Создаем приоритетную очередь (тип хранимого значения, в чем хранить, компаратор (по возрастанию)) для хранения вершин с минимальным расстоянием
        priority_queue<pair<double, int>, vector<pair<double, int>>, greater<pair<double, int>>> queue;
        // Добавляем стартовую вершину в очередь
        queue.push({ 0.0, startPosition });

        // Повторяем следующие шаги до тех пор, пока не будут посещены все вершины
        while (!queue.empty()) {
            // Извлекаем вершину с минимальным расстоянием из очереди
            double distance = queue.top().first;
            int position = queue.top().second;
            queue.pop();

            // Если мы достигли конечной вершины, возвращаем расстояния
            if (position == endPosition) {
                return distances;
            }

            // Для всех соседних вершин рассчитываем новое расстояние
            for (int j = 0; j < vertices.size(); ++j) {
                // Получаем вес ребра между текущей вершиной и соседней вершиной
                double weight = 0.0;
                if (position > j) {
                    weight = adjacencyMatrix[position][j];
                }
                else if (position < j) {
                    weight = adjacencyMatrix[j][position];
                }

                // Если ребро существует (т.е. вес не равен 0)
                if (weight != 0.0) {
                    // Рассчитываем новое расстояние до соседней вершины
                    double newDistance = distance + weight;

                    // Если новое расстояние меньше текущего расстояния, обновляем его
                    if (newDistance < distances[j]) {
                        distances[j] = newDistance;
                        // Добавляем соседнюю вершину в очередь
                        queue.push({ newDistance, j });
                    }
                }
            }
        }

        // Если мы не смогли найти путь до конечной вершины, бросаем исключение
        throw out_of_range("Путь между вершинами не найден.");
    }





    //// Процедура обхода графа в глубину
    //void depthFirstSearch(const T& startVertex) {
    //    // Находим позицию стартовой вершины в графе
    //    int startPosition = findVertex(startVertex);

    //    // Если вершина не найдена, бросаем исключение
    //    if (startPosition == -1) {
    //        throw out_of_range("Вершина не найдена в графе.");
    //    }

    //    // Создаем массив для отслеживания посещенных вершин
    //    vector<bool> visited(size, false);

    //    // Определяем функцию dfs, которая принимает позицию вершины
    //    function<void(int)> dfs = [&](int position) {
    //        // Помечаем вершину как посещенную
    //        visited[position] = true;
    //        // Выводим вершину на экран
    //        cout << vertices[position] << " ";

    //        // Для всех соседних вершин, которые еще не были посещены, вызываем функцию dfs рекурсивно
    //        for (int i = 0; i < size; ++i) {
    //            if (!visited[i] && adjacencyMatrix[position][i] != 0.0) {
    //                dfs(i);
    //            }
    //        }
    //        };

    //    // Вызываем функцию dfs для стартовой вершины
    //    dfs(startPosition);
    //    cout << endl;
    //}


    // Процедура обхода графа в глубину. Бросает исключение out_of_range если вершина не существует 
    void depthFirstSearch(const T& startVertex) const {
        // Находим позицию стартовой вершины в графе
        int startPosition = findVertex(startVertex);

        // Если вершина не найдена, бросаем исключение
        if (startPosition == -1) {
            throw out_of_range("Вершина не найдена в графе.");
        }

        // Создаем массив для отслеживания посещенных вершин
        vector<bool> visited(vertices.size(), false);

        // Определяем функцию dfs, которая принимает позицию вершины
        function<void(int)> dfs = [&](int position) {
            // Помечаем вершину как посещенную
            visited[position] = true;
            // Выводим вершину на экран
            cout << vertices[position] << " ";

            // Для всех соседних вершин, которые еще не были посещены, вызываем функцию dfs рекурсивно
            for (int i = 0; i < vertices.size(); ++i) {
                if (!visited[i]) {
                    double weight = 0.0;
                    if (position > i) {
                        weight = adjacencyMatrix[position][i];
                    }
                    else if (position < i) {
                        weight = adjacencyMatrix[i][position];
                    }
                    if (weight != 0.0) {
                        dfs(i);
                    }
                }
            }
            };
        // Вызываем функцию dfs для стартовой вершины
        dfs(startPosition);
        cout << endl;
    }

    //// Процедура обхода графа в ширину
    //void breadthFirstSearch(const T& startVertex) {
    //    // Находим позицию стартовой вершины в графе
    //    int startPosition = findVertex(startVertex);

    //    // Если вершина не найдена, бросаем исключение
    //    if (startPosition == -1) {
    //        throw out_of_range("Вершина не найдена в графе.");
    //    }

    //    // Создаем массив для отслеживания посещенных вершин и очередь для хранения вершин, которые нужно посетить
    //    vector<bool> visited(size, false);
    //    queue<int> queue;

    //    // Помечаем стартовую вершину как посещенную и добавляем ее в очередь
    //    queue.push(startPosition);
    //    visited[startPosition] = true;

    //    // Пока очередь не пуста, выполняем следующие действия
    //    while (!queue.empty()) {
    //        // Извлекаем вершину из очереди и выводим ее на экран
    //        int position = queue.front();
    //        queue.pop();
    //        cout << vertices[position] << " ";

    //        // Для всех соседних вершин, которые еще не были посещены, помечаем их как посещенные и добавляем их в очередь
    //        for (int i = 0; i < size; ++i) {
    //            if (!visited[i] && adjacencyMatrix[position][i] != 0.0) {
    //                queue.push(i);
    //                visited[i] = true;
    //            }
    //        }
    //    }
    //    cout << endl;
    //}

    // Процедура обхода графа в ширину. Бросает исключение out_of_range если вершина не существует 
    void breadthFirstSearch(const T& startVertex) const {
        // Находим позицию стартовой вершины в графе
        int startPosition = findVertex(startVertex);

        // Если вершина не найдена, бросаем исключение
        if (startPosition == -1) {
            throw out_of_range("Вершина не найдена в графе.");
        }

        // Создаем массив для отслеживания посещенных вершин и очередь для хранения вершин, которые нужно посетить
        vector<bool> visited(vertices.size(), false);
        queue<int> queue;

        // Помечаем стартовую вершину как посещенную и добавляем ее в очередь
        queue.push(startPosition);
        visited[startPosition] = true;

        // Пока очередь не пуста, выполняем следующие действия
        while (!queue.empty()) {
            // Извлекаем вершину из очереди и выводим ее на экран
            int position = queue.front();
            queue.pop();
            cout << vertices[position] << " ";

            // Для всех соседних вершин, которые еще не были посещены, помечаем их как посещенные и добавляем их в очередь
            for (int i = 0; i < vertices.size(); ++i) {
                if (!visited[i]) {
                    double weight = 0.0;
                    if (position > i) {
                        weight = adjacencyMatrix[position][i];
                    }
                    else if (position < i) {
                        weight = adjacencyMatrix[i][position];
                    }
                    if (weight != 0.0) {
                        queue.push(i);
                        visited[i] = true;
                    }
                }
            }
        }
        cout << endl;
    }


    //// Метод печати матрицы смежности
    //void printAdjacencyMatrix()
    //{
    //    cout << " " << setw(9);
    //    for (const T& name : vertices) {
    //        cout << name << setw(8) << " ";
    //    }
    //    cout << endl;
    //    for (size_t i = 0; i < adjacencyMatrix.size(); i++)
    //    {
    //        cout << vertices[i] << setw(8) << " ";
    //        for (size_t j = 0; j < adjacencyMatrix[0].size(); j++)
    //        {
    //            cout << adjacencyMatrix[i][j] << setw(8) << " ";
    //        }
    //        cout << endl;
    //    }
    //}

    // Метод печати матрицы смежности
    void printAdjacencyMatrix() const {
        cout << " " << setw(9);
        for (const T& name : vertices) {
            cout << name << setw(8) << " ";
        }
        cout << endl;
        for (size_t i = 0; i < adjacencyMatrix.size(); i++) {
            cout << vertices[i] << setw(8) << " ";
            for (size_t j = 0; j < i; j++) {
                cout << adjacencyMatrix[i][j] << setw(8) << " ";
            }
            for (size_t j = i; j < adjacencyMatrix[i].size(); j++) {
                cout << "0" << setw(8) << " ";
            }
            cout << endl;
        }
    }

    //// Метод обновления ребра между вершинами
    //void updateEdge(const T& startVertex, const T& endVertex, double weight) {
    //    int startPosition = findVertex(startVertex);
    //    int endPosition = findVertex(endVertex);

    //    if (startPosition == -1 || endPosition == -1) {
    //        throw out_of_range("Одна или обе вершины не найдены в графе.");
    //    }

    //    adjacencyMatrix[startPosition][endPosition] = weight;
    //    adjacencyMatrix[endPosition][startPosition] = weight;
    //}

    // Метод обновления ребра между вершинами. Бросает исключение out_of_range если вершина не существует 
    void updateEdge(const T& startVertex, const T& endVertex, double weight) {
        int startPosition = findVertex(startVertex);
        int endPosition = findVertex(endVertex);

        if (startPosition == -1 || endPosition == -1) {
            throw out_of_range("Одна или обе вершины не найдены в графе.");
        }

        if (startPosition > endPosition) {
            swap(startPosition, endPosition);
        }

        adjacencyMatrix[endPosition][startPosition] = weight;
    }


    //// Метод поиска всех путей между двумя вершинами
    //void findAllPaths(const T& startVertex, const T& endVertex) {
    //    // Находим позиции стартовой и конечной вершин в графе
    //    int startPosition = findVertex(startVertex);
    //    int endPosition = findVertex(endVertex);

    //    // Если одна или обе вершины не найдены, бросаем исключение
    //    if (startPosition == -1 || endPosition == -1) {
    //        throw out_of_range("Одна или обе вершины не найдены в графе.");
    //    }

    //    // Создаем массив для отслеживания посещенных вершин
    //    vector<bool> visited(size, false);
    //    // Создаем вектор для хранения пути
    //    vector<T> path;

    //    // Определяем функцию dfs, которая выполняет поиск в глубину
    //    function<void(int)> dfs = [&](int position) {
    //        // Помечаем вершину как посещенную
    //        visited[position] = true;
    //        // Добавляем вершину в путь
    //        path.push_back(vertices[position]);

    //        // Если вершина является конечной, выводим путь
    //        if (position == endPosition) {
    //            for (const auto& vertex : path) {
    //                cout << vertex << " ";
    //            }
    //            cout << endl;
    //        }
    //        else {
    //            // Для всех соседних вершин, которые еще не были посещены, вызываем функцию dfs
    //            for (int i = 0; i < size; ++i) {
    //                if (!visited[i] && adjacencyMatrix[position][i] != 0.0) {
    //                    dfs(i);
    //                }
    //            }
    //        }

    //        // Удаляем вершину из пути и помечаем ее как непосещенную
    //        path.pop_back();
    //        visited[position] = false;
    //        };

    //    // Вызываем функцию dfs для стартовой вершины
    //    dfs(startPosition);
    //}

    // Метод поиска всех путей между двумя вершинами. Бросает исключение out_of_range если вершина не существует 
    void findAllPaths(const T& startVertex, const T& endVertex) const {
        // Находим позиции стартовой и конечной вершин в графе
        int startPosition = findVertex(startVertex);
        int endPosition = findVertex(endVertex);

        // Если одна или обе вершины не найдены, бросаем исключение
        if (startPosition == -1 || endPosition == -1) {
            throw out_of_range("Одна или обе вершины не найдены в графе.");
        }

        // Создаем массив для отслеживания посещенных вершин
        vector<bool> visited(vertices.size(), false);
        // Создаем вектор для хранения пути
        vector<T> path;

        // Определяем функцию dfs, которая выполняет поиск в глубину
        function<void(int)> dfs = [&](int position) {
            // Помечаем вершину как посещенную
            visited[position] = true;
            // Добавляем вершину в путь
            path.push_back(vertices[position]);

            // Если вершина является конечной, выводим путь
            if (position == endPosition) {
                for (const auto& vertex : path) {
                    cout << vertex << " ";
                }
                cout << endl;
            }
            else {
                // Для всех соседних вершин, которые еще не были посещены, вызываем функцию dfs
                for (int i = 0; i < vertices.size(); ++i) {
                    if (!visited[i]) {
                        double weight = 0.0;
                        if (position > i) {
                            weight = adjacencyMatrix[position][i];
                        }
                        else if (position < i) {
                            weight = adjacencyMatrix[i][position];
                        }
                        if (weight != 0.0) {
                            dfs(i);
                        }
                    }
                }
            }

            // Удаляем вершину из пути и помечаем ее как непосещенную
            path.pop_back();
            visited[position] = false;
            };

        // Вызываем функцию dfs для стартовой вершины
        dfs(startPosition);
    }

    //// Метод вычисления степени вершины
    //int getVertexDegree(const T& vertex) {
    //    int position = findVertex(vertex);

    //    if (position == -1) {
    //        throw out_of_range("Вершина не найдена в графе.");
    //    }

    //    int degree = 0;
    //    for (int i = 0; i < size; ++i) {

    //        if (adjacencyMatrix[position][i] != 0.0) {
    //            degree++;
    //        }
    //        if (adjacencyMatrix[i][position] != 0.0) {
    //            degree++;
    //        }
    //    }

    //    return degree/2;
    //}

        // Метод вычисления степени вершины
    int getVertexDegree(const T& vertex) const {
        int position = findVertex(vertex);

        if (position == -1) {
            throw out_of_range("Вершина не найдена в графе.");
        }

        int degree = 0;
        for (int i = 0; i < vertices.size(); ++i) {
            if (i != position) {
                double weight = 0.0;
                if (position > i) {
                    weight = adjacencyMatrix[position][i];
                }
                else if (position < i) {
                    weight = adjacencyMatrix[i][position];
                }
                if (weight != 0.0) {
                    degree++;
                }
            }
        }

        return degree;
    }



    //// Метод вычисления суммы весов ребер
    //double getEdgesWeightSum() {
    //    double sum = 0.0;
    //    for (int i = 0; i < size; ++i) {
    //        for (int j = 0; j < size; ++j) {
    //            sum += adjacencyMatrix[i][j];
    //        }
    //    }

    //    return sum/2;
    //}

    // Метод вычисления суммы весов ребер
    double getEdgesWeightSum() const {
        double sum = 0.0;
        for (int i = 0; i < vertices.size(); ++i) {
            for (int j = 0; j < i; ++j) {
                sum += adjacencyMatrix[i][j];
            }
        }
        return sum;
    }

    // Процедура сохранения графа в формате GraphML
    void saveToGraphML(const string& filename) const {
        ofstream file(filename);
        if (!file.is_open()) {
            throw runtime_error("Не удалось открыть файл для записи.");
        }

        file << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << endl;
        file << "<graphml xmlns=\"http://graphml.graphdrawing.org/xmlns\"" << endl;
        file << "         xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"" << endl;
        file << "         xsi:schemaLocation=\"http://graphml.graphdrawing.org/xmlns" << endl;
        file << "         http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd\">" << endl;

        file << "  <key id=\"weight\" for=\"edge\" attr.name=\"weight\" attr.type=\"double\"\/>" << endl;
        file << "  <key id=\"label\" for=\"node\" attr.name=\"label\" attr.type=\"string\"\/>" << endl;
        // Описание графа
        file << "  <graph id=\"graph\" edgedefault=\"undirected\">" << endl;

        // Вершины
        for (size_t i = 0; i < vertices.size(); ++i) {
            file << "    <node id=\"" << i << "\">" << endl;
            //             file << "    <node id=\"" << vertices[i] << "\">" << endl;
            file << "      <data key=\"label\">" << vertices[i] << "</data>" << endl;
            file << "    </node>" << endl;
        }

        // Ребра
        for (size_t i = 0; i < vertices.size(); ++i) {
            for (size_t j = 0; j < i; ++j) {
                if (adjacencyMatrix[i][j] != 0.0) {
                    file << "    <edge id=\"" << i << "_" << j << "\" source=\"" << i << "\" target=\"" << j << "\" weight=\"" << adjacencyMatrix[i][j] << "\">" << endl;
                    file << "      <data key=\"weight\">" << adjacencyMatrix[i][j] << "</data>" << endl;
                    file << "    </edge>" << endl;
                }
            }
        }

        file << "  </graph>" << endl;
        file << "</graphml>" << endl;

        file.close();
    }



    //Тестирование класса
    static void runTests() {
        // Тестирование пустого графа
        WeightedGraphLegacy<string> emptyGraph;
        emptyGraph.printVertices();
        emptyGraph.printEdges();
        assert(emptyGraph.getEdgeCount() == 0);
        assert(emptyGraph.getVertexCount() == 0);
        try {
            emptyGraph.updateEdge("A", "B", 3.0);
            assert(false); // Если исключение не было отправлено, assert сработает
        }
        catch (const out_of_range& e) {
            // Если исключение было отправлено, assert не сработает
        }
        try {
            emptyGraph.findAllPaths("A", "B");
            assert(false); // Если исключение не было отправлено, assert сработает
        }
        catch (const out_of_range& e) {
            // Если исключение было отправлено, assert не сработает
        }
        try {
            emptyGraph.getVertexDegree("A") == 1;
            assert(false); // Если исключение не было отправлено, assert сработает
        }
        catch (const out_of_range& e) {
            // Если исключение было отправлено, assert не сработает
        }
        emptyGraph.getEdgesWeightSum() == 0;
        try {
            emptyGraph.findShortestPath("A", "B");
            assert(false); // Если исключение не было отправлено, assert сработает
        }
        catch (const out_of_range& e) {
            // Если исключение было отправлено, assert не сработает
        }

        try {
            emptyGraph.depthFirstSearch("A");
            assert(false); // Если исключение не было отправлено, assert сработает
        }
        catch (const out_of_range& e) {
            // Если исключение было отправлено, assert не сработает
        }

        try {
            emptyGraph.breadthFirstSearch("A");
            assert(false); // Если исключение не было отправлено, assert сработает
        }
        catch (const out_of_range& e) {
            // Если исключение было отправлено, assert не сработает
        }

        // Тестирование несвязного графа
        WeightedGraphLegacy<string> disconnectedGraph;
        disconnectedGraph.insertVertex("A");
        disconnectedGraph.insertVertex("B");
        disconnectedGraph.insertVertex("C");
        disconnectedGraph.insertVertex("D");

        assert(disconnectedGraph.getEdgeCount() == 0);
        assert(disconnectedGraph.getVertexCount() == 4);

        try {
            disconnectedGraph.findShortestPath("A", "C");
            assert(false); // Если исключение не было отправлено, assert сработает
        }
        catch (const out_of_range& e) {
            // Если исключение было отправлено, assert не сработает
        }

        disconnectedGraph.depthFirstSearch("A");
        disconnectedGraph.breadthFirstSearch("A");
        disconnectedGraph.findAllPaths("A", "C");

        assert(disconnectedGraph.getVertexDegree("A") == 0);
        assert(disconnectedGraph.getEdgesWeightSum() == 0);

        // Тестирование связного графа
        WeightedGraphLegacy<string> connectedGraph;
        connectedGraph.insertVertex("A");
        connectedGraph.insertVertex("B");
        connectedGraph.insertVertex("C");
        connectedGraph.insertVertex("D");
        connectedGraph.insertEdge("A", "B", 2.0);
        connectedGraph.insertEdge("B", "C", 3.0);
        connectedGraph.insertEdge("C", "D", 1.0);
        connectedGraph.insertEdge("D", "A", 4.0);

        assert(connectedGraph.getEdgeCount() == 4);
        assert(connectedGraph.getVertexCount() == 4);

        // Должен вызываться без исключений
        try {
            connectedGraph.findShortestPath("A", "C");
        }
        catch (const out_of_range& e) {
            assert(false); // Если исключение было отправлено, assert сработает
        }

        try {
            connectedGraph.depthFirstSearch("A");
        }
        catch (const out_of_range& e) {
            assert(false); // Если исключение было отправлено, assert сработает
        }

        try {
            connectedGraph.breadthFirstSearch("A");
        }
        catch (const out_of_range& e) {
            assert(false); // Если исключение было отправлено, assert сработает
        }

        connectedGraph.updateEdge("A", "B", 3.0);
        connectedGraph.findAllPaths("A", "C");
        assert(connectedGraph.getVertexDegree("A") == 2);
        assert(connectedGraph.getEdgesWeightSum() == 11.0);

        cout << "Все тесты пройдены успешно." << endl;
    }
};

WeightedGraphLegacy<string> createRandomGraph(size_t vertexCount, size_t maxEdgesPerVertex) {
    WeightedGraphLegacy<string> graph;

    // Генерируем названия вершин
    int nameLength = 1;
    while (pow(26, nameLength) < vertexCount) {
        nameLength++;
    }

    vector<string> vertexNames;
    for (int i = 0; i < vertexCount; i++) {
        string name;
        int num = i;
        for (int j = 0; j < nameLength; j++) {
            name = static_cast<char>('A' + num % 26) + name;
            num /= 26;
        }
        vertexNames.push_back(name);
    }

    // Добавляем вершины в граф
    for (const auto& name : vertexNames) {
        graph.insertVertex(name);
    }

    // Создаем случайные ребра
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> weightDistribution(1.0, 10.0);
    uniform_int_distribution<int> vertexDistribution(0, vertexCount - 1);

    for (int i = 0; i < vertexCount; i++) {
        int edgesCount = min(maxEdgesPerVertex, vertexCount - 1);
        for (int j = 0; j < edgesCount; j++) {
            int endVertexIndex;
            do {
                endVertexIndex = vertexDistribution(gen);
            } while (endVertexIndex == i);

            double weight = weightDistribution(gen);
            graph.insertEdge(vertexNames[i], vertexNames[endVertexIndex], weight);
        }
    }

    return graph;
}