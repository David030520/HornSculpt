// Included inside wideband_leaf.cu's internal namespace. These integrated
// signatures use exactly the canonical quadrature and scaled regular waves;
// no quadrature-point tableau or independent approximation is retained.
struct LeafFluxMomentSlot {
    unsigned node{}, face{};
};
struct LeafPanelMomentSlots {
    unsigned pressure[3]{absent, absent, absent};
    unsigned flux{absent};
};
__device__ std::size_t positiveMomentIndex(unsigned l, unsigned m) {
    return static_cast<std::size_t>(l) * (l + 1) / 2 + m;
}
__device__ Z momentOrder(const Z *moments, std::size_t base, unsigned l, int m) {
    const unsigned positive = m < 0 ? static_cast<unsigned>(-m) : static_cast<unsigned>(m);
    auto value = moments[base + positiveMomentIndex(l, positive)];
    // Real geometry, shape functions and k imply R_l,-m=(-1)^m conj(R_l,m),
    // including normal derivatives. The complex solution is never conjugated.
    if (m < 0)
        value = times(conjugate(value), positive % 2 ? -1 : 1);
    return value;
}
__global__ void buildLeafMoments(const LeafPanel *panels, const std::size_t *offsets,
                                 const LeafPanelMomentSlots *slots,
                                 const galerkin::TriangleQuadraturePoint *rule, unsigned ruleCount,
                                 unsigned order, double k, double radius, unsigned cachedPressure,
                                 unsigned cachedFlux, Z *normalMoments, Z *massMoments,
                                 Z *fluxMoments, unsigned first, std::size_t firstPanel,
                                 std::size_t panelCount, int *error) {
    extern __shared__ double radial[];
    const unsigned node = first + blockIdx.x;
    const auto count = static_cast<std::size_t>(order + 1) * (order + 2) / 2;
    // A P1 DOF and a (leaf, face) DP0 slot have one owning leaf. A physical
    // panel can occur in several leaves, with a different local origin in each.
    // Within a leaf, one m-column owns each coefficient of every output row.
    const auto panelsInLeaf = offsets[node + 1] - offsets[node];
    const auto relativeBegin = firstPanel < panelsInLeaf ? firstPanel : panelsInLeaf;
    const auto begin = offsets[node] + relativeBegin, remaining = panelsInLeaf - relativeBegin;
    const auto end = begin + (panelCount < remaining ? panelCount : remaining);
    // Consecutive slices retain each coefficient's complete panel/q sum order.
    // An exhausted leaf has an empty block-uniform interval, with no barriers.
    for (auto p = begin; p < end; ++p) {
        const auto entry = slots[p];
        const bool pressure = entry.pressure[0] < cachedPressure ||
                              entry.pressure[1] < cachedPressure ||
                              entry.pressure[2] < cachedPressure;
        if (!pressure && entry.flux >= cachedFlux)
            continue;
        const auto panel = panels[p];
        for (unsigned q = 0; q < ruleCount; ++q) {
            const auto quadrature = rule[q];
            const auto point = position(panel, quadrature);
            if (!threadIdx.x &&
                !regular::fillScaledRegularRadial(order, k * sqrt(dot(point, point)), radial))
                atomicExch(error, 1);
            __syncthreads();
            for (unsigned m = threadIdx.x; m <= order; m += blockDim.x) {
                regular::ScalarWave previous{}, current = regular::solidDiagonal(m, point, radius);
                for (unsigned l = m; l <= order; ++l) {
                    if (l > m) {
                        const auto next =
                            regular::solidNext(l, m, point, radius, previous, current);
                        previous = current;
                        current = next;
                    }
                    const auto wave =
                        regular::combineSolidRadial(current, l, point, k, radial[l], radial[l + 1]);
                    const auto coefficient = positiveMomentIndex(l, m);
                    if (pressure) {
                        const auto normal = derivative(wave, panel.normal);
                        for (unsigned local = 0; local < 3; ++local)
                            if (entry.pressure[local] < cachedPressure) {
                                const double weight = 2 * panel.area * quadrature.weight *
                                                      quadrature.barycentric[local];
                                const auto index = entry.pressure[local] * count + coefficient;
                                normalMoments[index] =
                                    add(normalMoments[index], times(normal, weight));
                                massMoments[index] =
                                    add(massMoments[index], times(wave.value, weight));
                            }
                    }
                    if (entry.flux < cachedFlux) {
                        const double weight = 2 * panel.area * quadrature.weight;
                        const auto index = entry.flux * count + coefficient;
                        fluxMoments[index] = add(fluxMoments[index], times(wave.value, weight));
                    }
                }
            }
            // In particular, m=order=128 is processed by thread zero's second
            // m-column before the shared radial table may be overwritten.
            __syncthreads();
        }
    }
}
template <bool pressure>
__global__ void
emitIntegratedMoments(const Z *moments, const LeafTest *incidences, const std::size_t *testOffsets,
                      const LeafFluxMomentSlot *fluxSlots, const std::size_t *nodeOffsets,
                      unsigned cached, unsigned order, unsigned storageOrder, const Z *input,
                      Z *output, std::size_t first, std::size_t total) {
    const auto within = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (within >= total)
        return;
    const auto flat = first + within, count = static_cast<std::size_t>(order + 1) * (order + 1),
               node = flat / count, coefficient = flat % count;
    const auto l = static_cast<unsigned>(sqrt(static_cast<double>(coefficient)));
    const int m =
        static_cast<int>(coefficient - static_cast<std::size_t>(l) * l) - static_cast<int>(l);
    const auto positiveCount = static_cast<std::size_t>(storageOrder + 1) * (storageOrder + 2) / 2;
    Z result{};
    for (auto slot = nodeOffsets[node]; slot < nodeOffsets[node + 1] && slot < cached; ++slot) {
        const unsigned column = pressure ? incidences[testOffsets[slot]].dof : fluxSlots[slot].face;
        const auto basis = momentOrder(moments, slot * positiveCount, l, m);
        result = add(result, multiply(input[column], conjugate(basis)));
    }
    output[flat] = add(output[flat], result);
}
__global__ void testIntegratedMoments(const Z *normalMoments, const Z *massMoments,
                                      const LeafTest *incidences, const std::size_t *offsets,
                                      unsigned order, unsigned storageOrder, double inverseK,
                                      const Z *input, Z *output, unsigned first, int *error) {
    __shared__ Z partial[threads];
    const unsigned target = first + blockIdx.x;
    const auto entry = incidences[offsets[target]];
    const auto count = static_cast<std::size_t>(order + 1) * (order + 1),
               positiveCount = static_cast<std::size_t>(storageOrder + 1) * (storageOrder + 2) / 2,
               base = static_cast<std::size_t>(target) * positiveCount;
    Z result{};
    for (auto coefficient = static_cast<std::size_t>(threadIdx.x); coefficient < count;
         coefficient += blockDim.x) {
        const auto l = static_cast<unsigned>(sqrt(static_cast<double>(coefficient)));
        const int m =
            static_cast<int>(coefficient - static_cast<std::size_t>(l) * l) - static_cast<int>(l);
        const auto d = momentOrder(normalMoments, base, l, m),
                   mass = momentOrder(massMoments, base, l, m);
        const auto tested =
            add(times(mass, -1), make_cuDoubleComplex(-d.y * inverseK, d.x * inverseK));
        result = add(
            result,
            multiply(tested, input[static_cast<std::size_t>(entry.node) * count + coefficient]));
    }
    partial[threadIdx.x] = result;
    __syncthreads();
    for (unsigned stride = blockDim.x / 2; stride; stride /= 2) {
        if (threadIdx.x < stride)
            partial[threadIdx.x] = add(partial[threadIdx.x], partial[threadIdx.x + stride]);
        __syncthreads();
    }
    if (!threadIdx.x) {
        const auto value = add(output[entry.dof], partial[0]);
        if (!isfinite(value.x) || !isfinite(value.y))
            atomicExch(error, 1);
        output[entry.dof] = value;
    }
}
