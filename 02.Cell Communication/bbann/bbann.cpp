#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <vector>
#include <set>
#include <omp.h>

struct Sorter
{
    size_t index;
    float value;
};

bool operator<(const Sorter& lhs, const Sorter& rhs)
{
    return lhs.value < rhs.value;
}

namespace py = pybind11;

py::tuple get_graph(py::array_t<float> affinity, std::vector<int> batch_key, int n_batches, int neighbors_within_batch)
{
    auto r = affinity.unchecked<2>();
    size_t I = r.shape(0);
    size_t n_ann = n_batches * neighbors_within_batch;
    std::vector<size_t> jj(I * n_ann);
    std::vector<float> val(I * n_ann);
    #pragma omp parallel for
    for(size_t i = 0; i < I; ++i)
    {
        std::vector<std::set<Sorter> *> batch_sorted(n_batches);
        for(size_t j = 0; j < n_batches; ++j)
            batch_sorted[j] = new std::set<Sorter>;
        std::set<Sorter> sorted;
        for (size_t j = 0; j < I; ++j)
        {
            struct Sorter ins = {j, r(i, j)};
            batch_sorted[batch_key[j]]->insert(ins);
        }
        for (size_t b = 0; b < n_batches; ++b)
        {
            size_t k = 0;
            for (std::set<Sorter>::reverse_iterator j = batch_sorted[b]->rbegin(); k < neighbors_within_batch; ++j, ++k)
            {
                struct Sorter ins = {j->index, j->value};
                sorted.insert(ins);
            }
        }
        size_t k = 0;
        for (std::set<Sorter>::reverse_iterator j = sorted.rbegin(); k < n_ann; ++j, ++k)
        {
            jj[i * n_ann + k] = j->index;
            val[i * n_ann + k] = j->value;
        }
        for (size_t j = 0; j < n_batches; ++j)
        {
            delete batch_sorted[j];
            batch_sorted[j] = NULL;
        }
    }
    py::array_t<int> indices = py::cast(jj);
    py::array_t<float> value = py::cast(val);
    return py::make_tuple(indices, value);
}

py::tuple calc_contrib_all(py::array_t<float> expression_a, py::array_t<float> expression_b)
{
    auto expr_a = expression_a.unchecked<2>();
    auto expr_b = expression_b.unchecked<2>();
    size_t na = expr_a.shape(0);
    size_t nb = expr_b.shape(0);
    size_t ng = expr_a.shape(1);
    std::vector<float> affivec(ng);
    std::vector<float> contrib(ng);
    #pragma omp parallel for
    for (size_t g = 0; g < ng; ++g)
        for (size_t i = 0; i < na; ++i)
            for (size_t j = 0; j < nb; ++j)
                affivec[g] += expr_a(i, g) * expr_b(j, g) / na / nb;
    float totalaffy = 0;
    for (size_t g = 0; g < ng; ++g)
        totalaffy += affivec[g];
    #pragma omp parallel for simd
    for (size_t g = 0; g < ng; ++g)
        contrib[g] = affivec[g] / totalaffy;
    return py::make_tuple(contrib, affivec);
}

py::tuple calc_contrib_ann(std::vector<int> row, std::vector<int> col, py::array_t<float> expression_a, py::array_t<float> expression_b)
{
    size_t n_pairs = row.size();
    auto expr_a = expression_a.unchecked<2>();
    auto expr_b = expression_b.unchecked<2>();
    size_t ng = expr_a.shape(1);
    std::vector<float> affinity(ng);
    std::vector<float> contrib(ng);
    #pragma omp parallel for
    for (size_t g = 0; g < ng; ++g)
        for (size_t m = 0; m < n_pairs; ++m)
            affinity[g] += expr_a(row[m], g) * expr_b(col[m], g) / n_pairs;
    float totalaffy = 0;
    for (size_t g = 0; g < ng; ++g)
        totalaffy += affinity[g];
    #pragma omp simd
    for (size_t g = 0; g < ng; ++g)
        contrib[g] = affinity[g] / totalaffy;
    return py::make_tuple(contrib, affinity);
}

py::array_t<int> connection(py::array_t<int> graph, std::vector<int> labels, int n_labels, std::vector<int> labels_n)
{
    auto ann = graph.unchecked<2>();
    py::array_t<double> connections({ n_labels, n_labels });
    auto conn = connections.mutable_unchecked<2>();
    size_t I = ann.shape(0);
    size_t J = ann.shape(1);
    for (size_t i = 0; i < I; ++i)
        for (size_t j = 0; j < J; ++j)
        {
            conn(labels[i], labels[ann(i, j)]) += 1;
            conn(labels[ann(i, j)], labels[i]) += 1;
        }
    for (size_t i = 0; i < n_labels; ++i)
        for (size_t j = 0; j < n_labels; ++j)
            conn(i, j) = conn(i, j) / labels_n[i] / labels_n[j];
    return connections;
            
}

PYBIND11_MODULE(bbann_internals, m) {
    m.def("get_graph", &get_graph, "Generate batch-balanced graph from affinity matrix.");
    m.def("calc_contrib_all", &calc_contrib_all, "Calculate contribution of gene pairs to affinity.");
    m.def("calc_contrib_ann", &calc_contrib_ann, "Calculate contribution of gene pairs to affinity.");
}
