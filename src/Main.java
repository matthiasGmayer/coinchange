import org.w3c.dom.ls.LSOutput;

import javax.imageio.ImageTypeSpecifier;
import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.concurrent.*;
import java.util.function.BiFunction;
import java.util.function.Function;
import java.util.function.Supplier;

public class Main {
    /*
    int[] m >= 0, distinct
    m[0]=1
    int k
    bestimme int[]n sodass sum m*n = k, sum n minimal
     */

    static boolean printSingle= false;

    public static void main(String[] args) {
//        normalTest(Main::alphaBeta,Main::recursiveGoodSpecialKgv,1,100,10000,true);

//        specificTest(new int[]{1,6,9,16,22,28,31,34,37,39,46,49,53,57,69,78,81,84,87,89,91,95,98,101,103,104,108,111}, 3468736);
//        lengthTest((m,k)->recursiveGoodSpecialLimit(m,k), 10000,1000);git
//        lengthTest((m,k)->recursiveGoodSpecialKgv(m,k), 10000,1000);
        lengthTest((m,k)->alphaBetaKGVs(m,k),10000,10000);
//        benchmarkTest((m,k)->alphaBetaKGVs(m,k,10),"benchHard.txt",true);
        speedMap.forEach((v,k)-> {
            int maxi=0;
            for (int i = 1; i < k.length; i++) {
                if(k[maxi]<k[i]) maxi=i;
            }
            System.out.println(v + " " +maxi);
        }
        );
//        System.out.println(t1);
//        System.out.println(t2);
//        benchmarkTest(Main::alphaBetaRaw,"bench2.txt",false);
    }
    private static long benchmarkTest(BiFunction<int[],Integer,int[]> func, String file, boolean testEquality){
        int bars = 20;
        for (int i = 0; i < bars; i++) {
            System.out.print("_");
        }
        System.out.println();
        var l = readBenchmark(file);
        int loops = l.size();
        long time = 0;
        int i=0;
        for (var data:
             l) {
            long t1 = System.nanoTime();
            var sol = func.apply(data.m,data.k);
            time += System.nanoTime()-t1;
            if(testEquality){
                int ourK = 0;
                for (int j = 0; j < sol.length; j++) {
                    ourK += sol[j];
                }
                if(ourK!=data.solN) System.out.println("Not Correct: "+(ourK- data.solN));
            }

            if(bars*i%loops==0)
                System.out.print(".");
            i++;
        }
        System.out.println();
        System.out.println("Took: "+time/1000000+" MilliSeconds");
        return time;
    }
    private static List<TestData> readBenchmark(String destin) {
        List<TestData> list = new ArrayList<>();
        File file = new File(destin);
        try {
            BufferedReader reader = new BufferedReader(new FileReader(file));
            String s1,s2,s3,s4;
            while((s1= reader.readLine())!=null) {
                s2 = reader.readLine();
                s3 = reader.readLine();
                s4 = reader.readLine();
                if(s2==null || s3 ==null|| s4 ==null) break;
                list.add(new TestData(s1,s2,s3,s4));
            }

        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
        return list;
    }
    private static class TestData{
        int k;
        int[]m;
        int solN;
        int[]sol;
        public TestData(String s1,String s2,String s3,String s4){
            k = Integer.parseInt(s1);
            m = new int[Integer.parseInt(s2)];
            var o =Arrays.stream(s3.split(" ")).map(s -> Integer.parseInt(s)).toArray();
            for (int i = 0; i < o.length; i++) {
                m[i]=(int)o[i];
            }
            var p= s4.split(" ");
            solN=Integer.parseInt(p[0]);
            sol = new int[m.length];
            for (int i = 0; i < m.length; i++) {
                sol[i] = Integer.parseInt(p[i+1]);
            }
            Arrays.sort(m);
        }
    }
    private static void specificTest(int[] ints, int i) {
        Arrays.sort(ints);
        Wrapper w = new Wrapper();
        var r = measure(()->recursiveGoodSpecialLimit(ints,i), w);
        print(r);
        System.out.println("Time: "+w.value);
    }

    public static int lengthTest(BiFunction<int[],Integer,int[]> func, int max, long timeLimit){
        int length = 2;
        while(true){
            int[] m = new int[length];
            m[0]=1;
            for (int j = 1; j < m.length; j++) {
                m[j]= (int)(Math.random()*max)+1;
            }
            int k = Arrays.stream(m).sum()+ (int)(Math.random()*100000);
            Arrays.sort(m);
            boolean changed = true;
            while(changed){
                changed = false;
                for (int j = 1; j < m.length; j++) {
                    if(m[j]<=m[j-1]){
                        m[j]++;
                        changed=true;
                    }
                }
            }

            final ExecutorService service = Executors.newSingleThreadExecutor();

            try {
                final Future<int[]> f = service.submit(() -> func.apply(m,k));
                f.get(timeLimit, TimeUnit.MILLISECONDS);
                System.out.print(length+"; ");
                if(length%30==0) System.out.println();
            } catch (final TimeoutException e) {
                System.out.println("Length "+length+" timed out");
                return length;
            } catch (final Exception e) {
                throw new RuntimeException(e);
            } finally {
                service.shutdown();
            }
            length++;
        }
    }
    public static void normalTest(BiFunction<int[],Integer,int[]> func,BiFunction<int[],Integer,int[]> func2, int loops, int coins, int max, boolean testEquality) {
        long sum1=0, sum2=0;
        int bars =20;
        for (int i = 0; i < bars; i++) {
            System.out.print("_");
        }
        System.out.println();
        int length = coins;
        for (int i = 0; i <loops; i++) {
//            int length = (int)(Math.random()*2)+3;

            int[] m = new int[length];
            m[0]=1;
            for (int j = 1; j < m.length; j++) {
                m[j]= (int)(Math.random()*max)+1;
            }

//            int[] m = new int[]{1,2,5,10,20,50,100,200,500};
//            int k = 78778695;
            int k = Arrays.stream(m).sum()+ (int)(Math.random()*10);
            Arrays.sort(m);
            boolean changed = true;
            while(changed){
                changed = false;
            for (int j = 1; j < m.length; j++) {
                if(m[j]<=m[j-1]){
                    m[j]++;
                    changed=true;
                }
            }
            }
            if(printSingle){
                System.out.print("M: ");
                print(m);
                System.out.println("K: "+k);
            }

            int[] finalM = m;
            Wrapper s1=new Wrapper(),s2=new Wrapper();
            var m1 = measure(()->func.apply(finalM,k), s1);
            var m2 =  measure(()->func2.apply(finalM,k), s2);
            if(testEquality&&(sumArray(m,m1)!=k||sumArray(m,m2)!=k||sumArray(m1)!=sumArray(m2))){
                System.out.print("M: ");
                print(m);
                System.out.println("K: "+k);
                print(m1);
                System.out.println(sumArray(m,m1)==k);
                print(m2);
                System.out.println(sumArray(m,m2)==k);
            }

            if(printSingle){
                print(m1);
                print(m2);
                System.out.println(s1.value+"\t\t\t"+s2.value);
                System.out.println((float)s1.value/s2.value +"\t\t\t"+(float)s2.value/s1.value);
            }else{
                if(bars*i%loops==0)
                    System.out.print(".");
            }
            sum1 += s1.value;
            sum2 += s2.value;
        }
        System.out.println();
        System.out.println("Summed");
        System.out.println(sum1+"\t\t\t"+sum2);
        System.out.println((float)sum1/sum2 +"\t\t\t"+(float)sum2/sum1);

    }
    public static void laufzeitTest(){
        long[] sum = new long[4];
        int loops=1000;
        int bars =20;
        for (int i = 0; i < bars; i++) {
            System.out.print("_");
        }
        System.out.println();
        for (int i = 0; i <loops; i++) {
//            int length = (int)(Math.random()*2)+3;
            int length = 5;
            int max = 100;

            int[] m = new int[length];
            m[0]=1;
            for (int j = 1; j < m.length; j++) {
                m[j]= (int)(Math.random()*max)+1;
            }

//            int[] m = new int[]{1,2,5,10,20,50,100,200,500};
//            int k = 78778695;
            int k = Arrays.stream(m).sum()+ (int)(Math.random()*100000000);
            Arrays.sort(m);
            boolean changed = true;
            while(changed){
                changed = false;
                for (int j = 1; j < m.length; j++) {
                    if(m[j]<=m[j-1]){
                        m[j]++;
                        changed=true;
                    }
                }
            }
            if(printSingle){
                System.out.print("M: ");
                print(m);
                System.out.println("K: "+k);
            }

            int[] finalM = m;
            Wrapper[] s=new Wrapper[4];
            for (int j = 0; j < s.length; j++) {
                s[j]=new Wrapper();
            }
            var m1 = recursiveGoodSpecialLimit(finalM,k,s);
                if(bars*i%loops==0)
                    System.out.print(".");

            for (int j = 0; j < sum.length; j++) {
                sum[j]+=s[j].value;
            }
        }
        System.out.println();
        System.out.println("Summed");
        for (int i = 0; i < sum.length; i++) {
            System.out.print(sum[i]/100000+" ");
        }
        System.out.println();
    }
    private static class Wrapper{
        long value;
    }
    public static<T> T measure(Supplier<T> s, Wrapper l){
        long t = System.nanoTime();
        var r = s.get();
        l.value= System.nanoTime()-t;
        return r;
    }
    public static int sumArray(int[]m, int[]n){
        int sum =0;
        for (int i = 0; i < m.length; i++) {
            sum += m[i] *n[i];
        }
        return sum;
    }
    public static int sumArray(int[]n){
        int sum= 0;
        for (int i = 0; i < n.length; i++) {
            sum += n[i];
        }
        return sum;
    }
    public static void print(int[] n) {
        int sum = 0;
        for (int i = 0; i < n.length; i++) {
            System.out.print(n[i] + " ");
            sum += n[i];
        }
        System.out.println("sum: " + sum);
    }
    public static int[] recursiveGoodSpecialLimit(int[] m, int k){
        if(m.length<5) return recursiveGoodSpecialKgv(m,k);

        int[] limit = new int[m.length-1];
        for (int i = 0; i < limit.length; i++) {
            limit[i]=kgv(m[i],m[i+1])/m[i]-1;
            for (int j = i+2; j < m.length; j++) {
                limit[i]=Math.min(limit[i],kgv(m[i],m[j])/m[i]-1);
            }
        }
        int loop = Math.min(m.length-1, m.length-1);
        int[][] a = new int[loop][];
        int maxK =0;
        for (int i = 0; i < loop; i++) {
            int[] array= new int[limit[i]+1];
            for (int j = 1; j <= limit[i] ; j++) {
                int val= m[i] * j;
                if(maxK < val) maxK= val;
                array[j] = val;
            }
            a[i]=array;
        }
        int[][] Am = dynamic(m,maxK+1);
        for (int i = 0; i < a.length; i++) {
            int lim = 0;
            for (int j = 0; j < a[i].length; j++) {
                lim=Math.max(lim,Am[a[i][j]][i]);
            }
            limit[i]=lim;
        }
        int[] maxDis = new int[m.length-1];
        maxDis[0]=m[1]-1;
        for (int i = 1; i < maxDis.length; i++) {
            maxDis[i]= maxDis[i-1]+limit[i]*m[i];
        }
        int[] best = new int[m.length];
        return recursiveGoodSpecialLimit(m,k,m.length-1,best, maxDis,limit);

    }
    public static int[] recursiveGoodSpecialLimit(int[] m, int k, Wrapper[] wrappers){

        var t1 = System.nanoTime();
        int[] limit = new int[m.length-1];
        for (int i = 0; i < limit.length; i++) {
            limit[i]=kgv(m[i],m[i+1])/m[i]-1;
            for (int j = i+2; j < m.length; j++) {
                limit[i]=Math.min(limit[i],kgv(m[i],m[j])/m[i]-1);
            }
        }
        int[] temp = Arrays.copyOf(limit, limit.length);
        int loop = Math.min(m.length-1, m.length-1);
        int[][] a = new int[loop][];
        int maxK =0;
        for (int i = 0; i < loop; i++) {
//            int lim=0;
//            for (int j = 1; j <= limit[i] ; j++) {
//                lim=Math.max(lim,recursiveGoodSpecialKgv(m,m[i]*j)[i]);
//            }
//            limit[i]=lim;
            int[] array= new int[limit[i]+1];
            for (int j = 1; j <= limit[i] ; j++) {
                int val= m[i] * j;
                if(maxK < val) maxK= val;
                array[j] = val;
            }
            a[i]=array;
        }
        wrappers[0].value = System.nanoTime()-t1;
        //Dynamisch
        t1 = System.nanoTime();
        int[][] Am = dynamic(m,maxK+1);
        wrappers[1].value = System.nanoTime()-t1;
        t1 = System.nanoTime();

        for (int i = 0; i < a.length; i++) {
            int lim = 0;
            for (int j = 0; j < a[i].length; j++) {
                lim=Math.max(lim,Am[a[i][j]][i]);
            }
            limit[i]=lim;
        }

        int[] maxDis = new int[m.length-1];
        maxDis[0]=m[1]-1;
        for (int i = 1; i < maxDis.length; i++) {
            maxDis[i]= maxDis[i-1]+limit[i]*m[i];
        }
        int[] best = new int[m.length];
        wrappers[2].value = System.nanoTime()-t1;

        t1 = System.nanoTime();

        var ret = recursiveGoodSpecialLimit(m,k,m.length-1,best, maxDis,limit);
        wrappers[3].value = System.nanoTime()-t1;

        return ret;
    }
    public static int[] dynamicSingle(int[]m,int k){
        return dynamic(m,k)[k];
    }
    private static int[][] dynamic(int[] m, int k) {
        int[] minSum= new int[k+1];
        int[][] minM = new int[k+1][m.length];
        for (int i = 1; i < minM.length; i++) {
            boolean isM=false;
            for (int j = 0; j < m.length; j++) {
                if(i==m[j]){
                    minSum[i]=1;
                    minM[i][j]=1;
                    isM = true;
                    break;
                }
            }
            if(!isM){
                int min = Integer.MAX_VALUE;
                int[] mM = new int[m.length];
                for (int j = 0; j<m.length && m[j] < i; j++) {
                    int sum = minSum[m[j]]+minSum[i-m[j]];
                    if(sum < min){
                        min= sum;
                        for (int l = 0; l < mM.length; l++) {
                            mM[l]=minM[m[j]][l]+minM[i-m[j]][l];
                        }
                    }
                }
                minSum[i]=min;
                minM[i] =mM;
            }
        }
        return minM;
    }

    public static int[] recursiveGoodSpecialLimit(int[] m, int k, Wrapper wrapper){
        int[] limit = new int[m.length-1];
        for (int i = 0; i < limit.length; i++) {
            limit[i]=kgv(m[i],m[i+1])/m[i]-1;
            for (int j = i+2; j < m.length; j++) {
                limit[i]=Math.min(limit[i],kgv(m[i],m[j])/m[i]-1);
            }
        }
        int[] temp = Arrays.copyOf(limit, limit.length);
        int loop = Math.min(m.length-1, m.length-1);
        for (int i = 1; i < loop; i++) {
            int lim=0;
            for (int j = 1; j <= limit[i] ; j++) {
                lim=Math.max(lim,recursiveGoodSpecialKgv(m,m[i]*j)[i]);
            }
            limit[i]=lim;
        }
        int[] maxDis = new int[m.length-1];
        maxDis[0]=m[1]-1;
        for (int i = 1; i < maxDis.length; i++) {
            maxDis[i]= maxDis[i-1]+limit[i]*m[i];
        }
        int[] best = new int[m.length];
        return measure(()->recursiveGoodSpecialLimit(m,k,m.length-1,best, maxDis,limit),wrapper);
    }

    public static int[] recursiveGoodSpecialLimit(int[] m, int k, int level, int[] best, int[] maxDis, int[] limit){
        int[] n = best;
        int min = Integer.MAX_VALUE;
        if(level ==0){
            best[0] = k;
        }else{
            int start = (Math.max(k-maxDis[level-1],0)+m[level]-1)/m[level];
            int end = k/m[level];
            if(level<m.length-1){
                end= Math.min(end,limit[level]);
            }
            for (int i = end; i >= start; i--) {
                best = recursiveGoodSpecialLimit(m,k-i*m[level], level-1, best, maxDis,limit);
                int sum=i;
                for (int j = 0; j < level; j++) {
                    sum += best[j];
                }
                if(sum<min){
                    best[level]=i;
                    min = sum;
                    n= best;
                    best = new int[m.length];
                }
            }
        }
        return n;
    }
    public static int[] recursiveGoodSpecialKgv(int[] m, int k){
        int[] kgv = new int[m.length-1];
        for (int i = 0; i < kgv.length; i++) {
            kgv[i]=kgv(m[i],m[i+1])/m[i]-1;
            for (int j = i+2; j < m.length; j++) {
                kgv[i]=Math.min(kgv[i],kgv(m[i],m[j])/m[i]-1);
            }
        }
        int[] maxDis = new int[m.length-1];
        maxDis[0]=m[1]-1;
        for (int i = 1; i < maxDis.length; i++) {
            maxDis[i]= maxDis[i-1]+kgv[i]*m[i];
        }
        int[] best = new int[m.length];
        return recursiveGoodSpecialKgv(m,k,m.length-1,best, maxDis,kgv);
    }
    public static int[] recursiveGoodSpecialKgv(int[] m, int k, int level, int[] best, int[] maxDis, int[] kgv){
        int[] n = best;
        int min = Integer.MAX_VALUE;
        if(level ==0){
            best[0] = k;
        }else{
            int start = (Math.max(k-maxDis[level-1],0)+m[level]-1)/m[level];
            int end = k/m[level];
            if(level<m.length-1){
                end= Math.min(end,kgv[level]);
            }
            for (int i = end; i >= start; i--) {
                best = recursiveGoodSpecialKgv(m,k-i*m[level], level-1, best, maxDis,kgv);
                int sum=i;
                for (int j = 0; j < level; j++) {
                    sum += best[j];
                }
                if(sum<min){
                    best[level]=i;
                    min = sum;
                    n= best;
                    best = new int[m.length];
                }
            }
        }
        return n;
    }
    public static int kgv(int z1, int z2){
        int z1temp = z1;
        int z2temp = z2;
        while (z1temp != z2temp) {
            if (z1temp < z2temp) {
                z1temp += z1;
            } else {
                z2temp+= z2;
            }
        }
        return z1temp;
    }
    public static int[] alphaBetaLimit(int[] m, int k){
        int[] limit = new int[m.length-1];
        for (int i = 0; i < limit.length; i++) {
            limit[i]=kgv(m[i],m[i+1])/m[i]-1;
            for (int j = i+2; j < m.length; j++) {
                limit[i]=Math.min(limit[i],kgv(m[i],m[j])/m[i]-1);
            }
        }
        int loop = Math.min(m.length-1, m.length-1);
        int[][] a = new int[loop][];
        int maxK =0;
        for (int i = 0; i < loop; i++) {
            int[] array= new int[limit[i]+1];
            for (int j = 1; j <= limit[i] ; j++) {
                int val= m[i] * j;
                if(maxK < val) maxK= val;
                array[j] = val;
            }
            a[i]=array;
        }
        int[][] Am = dynamic(m,maxK+1);
        for (int i = 0; i < a.length; i++) {
            int lim = 0;
            for (int j = 0; j < a[i].length; j++) {
                lim=Math.max(lim,Am[a[i][j]][i]);
            }
            limit[i]=lim;
        }
        int[] maxDis = new int[m.length-1];
        maxDis[0]=m[1]-1;
        for (int i = 1; i < maxDis.length; i++) {
            maxDis[i]= maxDis[i-1]+limit[i]*m[i];
        }
        return alphaBeta(m,k,m.length-1, maxDis,limit,new int[]{Integer.MAX_VALUE},0);
    }
    static long t1=0;
    static long t2=0;
    public static int[] alphaBetaKGVs(int[]m,int k){
        long t = System.nanoTime();
        int[] kgv = new int[m.length];
        for (int i = 0; i < kgv.length-1; i++) {
            kgv[i]=kgv(m[i],m[i+1])/m[i]-1;
            for (int j = i+2; j < m.length; j++) {
                kgv[i]=Math.min(kgv[i],kgv(m[i],m[j])/m[i]-1);
            }
        }
//        for (int i = 0; i < kgv.length-1; i++) {
//            kgv[i]=kgv(m[i],m[i+1])/m[i]-1;
//        }
        kgv[kgv.length-1]=k/m[m.length-1];
        int[] maxDis = new int[m.length-1];
        maxDis[0]=m[1]-1;
        for (int i = 1; i < maxDis.length; i++) {
            maxDis[i]= maxDis[i-1]+kgv[i]*m[i];
        }
        t1 += System.nanoTime()-t;
        t=System.nanoTime();
        var r =  alphaBeta(m,k,m.length-1, maxDis,kgv, new int[]{Integer.MAX_VALUE},0);
        t2 += System.nanoTime()-t;
        return r;
    }

    public static int getBeta(int[] m, int k,int level,int[] limit,int used,int[]alpha){
        int beta = used;
        for(int i=level; i> 0; i--){
            int max = (int)Math.ceil((double)(k)/(double)(m[i]));
            max = Math.min(max,limit[i]);
            beta += max;
            k -= max * m[i];
            if(k <=0|| alpha[0] <= beta){return beta;}
        }
        beta += k;
        return beta;
    }
    public static int[] alphaBetaKGVs(int[]m,int k,int howMany){
        long t = System.nanoTime();
        int[] kgv = new int[m.length];
        for (int i = 0; i < kgv.length-1; i++) {
            kgv[i]=kgv(m[i],m[i+1])/m[i]-1;
            for (int j = i+2; j < Math.min(m.length,howMany+i+2); j++) {
                kgv[i]=Math.min(kgv[i],kgv(m[i],m[j])/m[i]-1);
            }
        }
//        for (int i = 0; i < kgv.length-1; i++) {
//            kgv[i]=kgv(m[i],m[i+1])/m[i]-1;
//        }
        kgv[kgv.length-1]=k/m[m.length-1];
        int[] maxDis = new int[m.length-1];
        maxDis[0]=m[1]-1;
        for (int i = 1; i < maxDis.length; i++) {
            maxDis[i]= maxDis[i-1]+kgv[i]*m[i];
        }
        t1 += System.nanoTime()-t;
        t=System.nanoTime();
        var r =  alphaBeta(m,k,m.length-1, maxDis,kgv, new int[]{Integer.MAX_VALUE},0);
        t2 += System.nanoTime()-t;
        return r;
    }
    public static int[] alphaBetaKGV(int[]m,int k){
        long t = System.nanoTime();
        int[] kgv = new int[m.length];
        for (int i = 0; i < kgv.length-1; i++) {
            kgv[i]=kgv(m[i],m[i+1])/m[i]-1;
        }
        kgv[kgv.length-1]=k/m[m.length-1];
        int[] maxDis = new int[m.length-1];
        maxDis[0]=m[1]-1;
        for (int i = 1; i < maxDis.length; i++) {
            maxDis[i]= maxDis[i-1]+kgv[i]*m[i];
        }
        t1 += System.nanoTime()-t;
        t=System.nanoTime();
        var r =  alphaBeta(m,k,m.length-1, maxDis,kgv, new int[]{Integer.MAX_VALUE},0);
        t2 += System.nanoTime()-t;
        return r;
    }
    public static int[] alphaBetaNeighbour(int[]m,int k){
        long t = System.nanoTime();
        int[] kgv = new int[m.length];
        for (int i = 0; i < kgv.length-1; i++) {
            kgv[i]=m[i]*m[i+1]/m[i]-1;
        }
        kgv[kgv.length-1]=k/m[m.length-1];
        int[] maxDis = new int[m.length-1];
        maxDis[0]=m[1]-1;
        for (int i = 1; i < maxDis.length; i++) {
            maxDis[i]= maxDis[i-1]+kgv[i]*m[i];
        }
        t1 += System.nanoTime()-t;
        t=System.nanoTime();
        var r =  alphaBeta(m,k,m.length-1, maxDis,kgv, new int[]{Integer.MAX_VALUE},0);
        t2 += System.nanoTime()-t;
        return r;
    }
    public static int[] alphaBeta(int[] m, int k, int level,int[] maxDis,int[] limit,int []alpha, int used){
        int[] sol = new int[level+1];
        //Wenn die münze k teil haben wir eine Lösung
        int div = k/m[level];
        if(div*m[level]==k) {
            if(used+div<alpha[0]) {
//                System.out.println(alpha[0]+ "  " +used + div);
                alpha[0] = used + div;
            }
            sol[level]=div;
            return sol;
        }
        int min = Integer.MAX_VALUE;
        int minIndex = (Math.max(k-maxDis[level-1],0)+m[level]-1)/m[level];
        int maxIndex = Math.min(k/m[level],limit[level]);
        int beta = used+(int)Math.ceil((double)k/(double)m[level]);
//        int beta= getBeta(m,k,level,limit,used,alpha);
//        if(beta!=beta2)
//            System.out.println(beta-beta2);
         if(alpha[0]<=beta){
            sol[0]=k;
            return sol;
        }
        //alle möglichkeiten durchprobieren
        for (int i = maxIndex; i >= minIndex; i--) {
            int currentValue = i*m[level];
            int[] partielSol = alphaBeta(m,k-currentValue,level-1, maxDis, limit, alpha,used+i);
            int sum = i;
            for (int j = 0; j < partielSol.length; j++) {
                sum+=partielSol[j];
            }
            if(sum<min){
                min=sum;
                sol[level]=i;
                for (int j = 0; j < partielSol.length; j++) {
                    sol[j]=partielSol[j];
                }
            }
            if(alpha[0]<=beta){
                break;

            }
        }
        return sol;
    }
    public static int[] alphaBetaRaw(int[] m, int k){
        return alphaBetaRaw(m,k,m.length-1,new int[]{Integer.MAX_VALUE},0);
    }
    public static int[] alphaBetaRaw(int[] m, int k, int level,int []alpha, int used){
        int[] sol = new int[level+1];
        //Wenn die münze k teil haben wir eine Lösung
        int div = k/m[level];
        if(div*m[level]==k) {
            if(used+div<alpha[0]) {
                alpha[0] = used + div;
            }
            sol[level]=div;
            return sol;
        }
        int min = Integer.MAX_VALUE;
        int minIndex = 0;
        int maxIndex = k/m[level];
        int beta = used+(int)Math.ceil((double)k/(double)m[level]);
        if(alpha[0]<=beta){
            sol[0]=k;
            return sol;
        }
        //alle möglichkeiten durchprobieren
        for (int i = maxIndex; i >= minIndex; i--) {
            int currentValue = i*m[level];
            int[] partielSol = alphaBetaRaw(m,k-currentValue,level-1, alpha,used+i);
            int sum = i;
            for (int j = 0; j < partielSol.length; j++) {
                sum+=partielSol[j];
            }
            if(sum<min){
                min=sum;
                sol[level]=i;
                for (int j = 0; j < partielSol.length; j++) {
                    sol[j]=partielSol[j];
                }
            }
            if(alpha[0]<=beta){
                break;

            }
        }
        return sol;
    }
    public static int[] alphaBetaLimitAlpha(int[]m,int k){
        long t = System.nanoTime();
        int[] limit = new int[m.length];
        for (int i = 0; i < limit.length-1; i++) {
            limit[i]=kgv(m[i],m[i+1])/m[i]-1;
            for (int j = i+2; j < m.length; j++) {
                limit[i]=Math.min(limit[i],kgv(m[i],m[j])/m[i]-1);
            }
        }
        limit[limit.length-1]=k/m[m.length-1];
        int[] maxDis = new int[m.length-1];
        maxDis[0]=m[1]-1;
        for (int i = 1; i < maxDis.length; i++) {
            maxDis[i]= maxDis[i-1]+limit[i]*m[i];
        }
        int maxK =0;
        int[][] a = new int[m.length-1][];
        for (int i = 0; i < m.length-1; i++) {
            int[] array= new int[limit[i]+1];
            for (int j = 1; j <= limit[i] ; j++) {
                int val= m[i] * j;
                if(maxK < val) maxK= val;
                array[j] = val;
            }
            a[i]=array;
        }
        HashMap<Integer,int[]> map = new HashMap<>();
        for (int i = 0; i < a.length; i++) {
            for (int j = 0; j < a[i].length; j++) {
                map.put(a[i][j],alphaBeta(m,a[i][j],m.length-1,maxDis,limit,new int[]{Integer.MAX_VALUE},0));
            }
        }
        for (int i = 0; i < a.length; i++) {
            int lim = 0;
            for (int j = 0; j < a[i].length; j++) {
                lim=Math.max(lim,map.get(a[i][j])[i]);
            }
            limit[i]=lim;
        }

        for (int i = 1; i < maxDis.length; i++) {
            maxDis[i]= maxDis[i-1]+limit[i]*m[i];
        }
        t1 += System.nanoTime()-t;
        t=System.nanoTime();
        var r=  alphaBeta(m,k,m.length-1, maxDis,limit, new int[]{Integer.MAX_VALUE},0);
        t2 += System.nanoTime()-t;
        return r;
    }
    static HashMap<Integer, int[]> speedMap = new HashMap<>();
    public static int[] alphaBetaCombined(int[]m,int k){
        long t1,t2,t3;
        t1=System.nanoTime();
        var r = alphaBetaNeighbour(m, k);
        t1 = System.nanoTime()-t1;
        t2= System.nanoTime();
        var r2 = alphaBetaKGV(m,k);
        t2 = System.nanoTime()-t2;
        t3=System.nanoTime();
        var r3=alphaBetaKGVs(m,k);
        t3 = System.nanoTime()-t3;
        var ma = speedMap.get(m.length);
        if(ma == null) speedMap.put(m.length,ma = new int[3]);
        if(t1>=t2&&t1>=t3) ma[0]++;
        else if(t2>=t1&&t2>=t3) ma[1]++;
        else ma[2]++;
//        int l = m.length;
//        if(l<30) return alphaBetaNeighbour(m,k);
//        if(l<100) return alphaBetaKGV(m,k);
//        return alphaBetaKGVs(m,k);
        return r;
    }
}
