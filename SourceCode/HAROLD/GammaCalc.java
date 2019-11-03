package cluster_rg;

import org.apache.commons.math3.special.Gamma;
import org.cache2k.Cache;
import org.cache2k.Cache2kBuilder;
import org.cache2k.integration.CacheLoader;

public abstract class GammaCalc {

    public abstract double logGamma(final double x);

    public static GammaCalc get(final int cacheSize) {
        if (cacheSize == 0) {
            return new GammaCalc() {
                @Override
                public final double logGamma(final double x) {
                    return Gamma.logGamma(x);
                }
            };
        } else {
            return new GammaCalc() {
                private final Cache<Double, Double> logGammaCache = new Cache2kBuilder<Double, Double>() { }
                        .name("logGamma")
                        .entryCapacity(cacheSize)
                        .loader(new CacheLoader<Double, Double>() {
                            @Override
                            public Double load(final Double x) {
                                return Gamma.logGamma(x);
                            }
                        })
                        .build();

                @Override
                public final double logGamma(final double x) {
                    return this.logGammaCache.get(x);
                }
            };
        }
    }
}

