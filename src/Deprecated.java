import java.util.Arrays;
import java.util.HashMap;

public class Deprecated {
    public static int[] recursiveGoodSpecialLimitHash(int[] m, int k){

        int[] limit = new int[m.length-1];
        for (int i = 0; i < limit.length; i++) {
            limit[i]=Main.kgv(m[i],m[i+1])/m[i]-1;
            for (int j = i+2; j < m.length; j++) {
                limit[i]=Math.min(limit[i],Main.kgv(m[i],m[j])/m[i]-1);
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
        int[][] Am = dynamicHash(m,maxK+1,a);
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
        return Main.recursiveGoodSpecialLimit(m,k,m.length-1,best, maxDis,limit);

    }
    private static int[][] dynamicHash(int[] m, int k, int[][] a){
        HashMap<Integer, int[]> map = new HashMap<>();
        HashMap<Integer, Integer> sumMap = new HashMap<>();
        int[][] minM = new int[k+1][m.length];
        for (int i = 0; i < a.length; i++) {
            for (int j = 0; j < a[i].length; j++) {
                minM[a[i][j]] = dynamicRecursiveHash(m,a[i][j],map,sumMap);
            }
        }
        for (int i = 0; i < k; i++) {
            if(!sumMap.keySet().contains(i)) System.out.println(i);
        }
        return minM;
    }
    private static int[] dynamicRecursiveHash(int[] m, int k, HashMap<Integer, int[]> map, HashMap<Integer, Integer> sumMap){
        var r = map.get(k);
        if(r!=null) return r;
        for (int j = 0; j < m.length; j++) {
            if(k==m[j]){
                sumMap.put(k,1);
                int[] min= new int[m.length];
                min[j]=1;
                map.put(k,min);
                return min;
            }
        }
        int min = Integer.MAX_VALUE;
        int[] mM = new int[m.length];
        for (int j = 0; j<m.length && m[j] < k; j++) {
            int sum = dynamicRecursiveHashSum(m,m[j],map,sumMap)+dynamicRecursiveHashSum(m,k-m[j],map,sumMap);
            if(sum < min){
                min= sum;
                for (int l = 0; l < mM.length; l++) {
                    mM[l]=dynamicRecursiveHash(m,m[j],map,sumMap)[l]+dynamicRecursiveHash(m,k-m[j],map,sumMap)[l];
                }
            }
        }
        sumMap.put(k,min);
        map.put(k,mM);
        return mM;
    }
    private static int dynamicRecursiveHashSum(int[] m, int k, HashMap<Integer, int[]> map, HashMap<Integer, Integer> sumMap){
        var r = sumMap.get(k);
        if(r!=null) return r;
        for (int j = 0; j < m.length; j++) {
            if(k==m[j]){
                sumMap.put(k,1);
                int[] min= new int[m.length];
                min[j]=1;
                map.put(k,min);
                return 1;
            }
        }
        int min = Integer.MAX_VALUE;
        int[] mM = new int[m.length];
        for (int j = 0; j<m.length && m[j] < k; j++) {
            int sum = dynamicRecursiveHashSum(m,m[j],map,sumMap)+dynamicRecursiveHashSum(m,k-m[j],map,sumMap);
            if(sum < min){
                min= sum;
                for (int l = 0; l < mM.length; l++) {
                    mM[l]=dynamicRecursiveHash(m,m[j],map,sumMap)[l]+dynamicRecursiveHash(m,k-m[j],map,sumMap)[l];
                }
            }
        }
        sumMap.put(k,min);
        map.put(k,mM);
        return min;
    }
    public static int[] recursiveGoodSpecialLimitHashHalfed(int[] m, int k){

        int[] limit = new int[m.length-1];
        for (int i = 0; i < limit.length; i++) {
            limit[i]=Main.kgv(m[i],m[i+1])/m[i]-1;
            for (int j = i+2; j < m.length; j++) {
                limit[i]=Math.min(limit[i],Main.kgv(m[i],m[j])/m[i]-1);
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
        int[] maxDis = new int[m.length-1];
        maxDis[0]=m[1]-1;
        for (int i = 1; i < maxDis.length; i++) {
            maxDis[i]= maxDis[i-1]+limit[i]*m[i];
        }
        int[][] Am = dynamicHash(m,maxK+1,a, maxDis);

        for (int i = 0; i < a.length; i++) {
            int lim = 0;
            for (int j = 0; j < a[i].length; j++) {
                lim=Math.max(lim,Am[a[i][j]][i]);
            }
            limit[i]=lim;
        }
        maxDis = new int[m.length-1];
        maxDis[0]=m[1]-1;
        for (int i = 1; i < maxDis.length; i++) {
            maxDis[i]= maxDis[i-1]+limit[i]*m[i];
        }
        int[] best = new int[m.length];
        return Main.recursiveGoodSpecialLimit(m,k,m.length-1,best, maxDis,limit);

    }
    private static int[][] dynamicHash(int[] m, int k, int[][] a, int [] maxDis){
        HashMap<Integer, int[]> map = new HashMap<>();
        HashMap<Integer, Integer> sumMap = new HashMap<>();
        int[][] minM = new int[k+1][m.length];
        for (int i = 0; i < a.length; i++) {
            for (int j = 0; j < a[i].length; j++) {
                minM[a[i][j]] = dynamicRecursiveHash(m,a[i][j],map,sumMap, maxDis);
            }
        }
        for (int i = 0; i < k; i++) {
            if(!sumMap.keySet().contains(i)) System.out.println(i);
        }
        return minM;
    }
    private static int[] dynamicRecursiveHash(int[] m, int k, HashMap<Integer, int[]> map, HashMap<Integer, Integer> sumMap, int[]maxDis){
        var r = map.get(k);
        if(r!=null) return r;
        for (int j = 0; j < m.length; j++) {
            if(k==m[j]){
                sumMap.put(k,1);
                int[] min= new int[m.length];
                min[j]=1;
                map.put(k,min);
                return min;
            }
        }
        int min = Integer.MAX_VALUE;
        int[] mM = new int[m.length];
        int minint = (Math.max(k-maxDis[maxDis.length-1],0)+m[maxDis.length]-1)/m[maxDis.length];
        for (int j = 0; j<m.length && m[j] < k; j++) {
            int sum = dynamicRecursiveHashSum(m,m[j]+minint,map,sumMap,maxDis)+dynamicRecursiveHashSum(m,k-m[j]-minint,map,sumMap,maxDis);
            if(sum < min){
                min= sum;
                for (int l = 0; l < mM.length; l++) {
                    mM[l]=dynamicRecursiveHash(m,m[j]+minint,map,sumMap,maxDis)[l]+dynamicRecursiveHash(m,k-m[j]-minint,map,sumMap,maxDis)[l];
                }
            }
        }
        sumMap.put(k,min);
        map.put(k,mM);
        return mM;
    }
    private static int dynamicRecursiveHashSum(int[] m, int k, HashMap<Integer, int[]> map, HashMap<Integer, Integer> sumMap, int[] maxDis){
        var r = sumMap.get(k);
        if(r!=null) return r;
        for (int j = 0; j < m.length; j++) {
            if(k==m[j]){
                sumMap.put(k,1);
                int[] min= new int[m.length];
                min[j]=1;
                map.put(k,min);
                return 1;
            }
        }
        int min = Integer.MAX_VALUE;
        int[] mM = new int[m.length];
        int minint = (Math.max(k-maxDis[maxDis.length-1],0)+m[maxDis.length]-1)/m[maxDis.length];
        for (int j = 0; j<m.length && m[j] < k; j++) {
            int sum = dynamicRecursiveHashSum(m,m[j]+minint,map,sumMap,maxDis)+dynamicRecursiveHashSum(m,k-m[j]-minint,map,sumMap,maxDis);
            if(sum < min){
                min= sum;
                for (int l = 0; l < mM.length; l++) {
                    mM[l]=dynamicRecursiveHash(m,m[j]+minint,map,sumMap,maxDis)[l]+dynamicRecursiveHash(m,k-m[j]-minint,map,sumMap,maxDis)[l];
                }
            }
        }
        sumMap.put(k,min);
        map.put(k,mM);
        return min;
    }
    public static int[] recursiveGoodSpecial2(int[] m, int k){
        int[] m2sum = new int[m.length];
        for (int i = 1; i < m2sum.length; i++) {
            m2sum[i]=m2sum[i-1]+m[i-1]*m[i];
        }
        int[] maxDis = new int[m.length-1];
        maxDis[0]=m[1]-1;
        for (int i = 1; i < maxDis.length; i++) {
            maxDis[i]= maxDis[i-1]+(m[i+1]-1)*m[i];
        }
        int[] best = new int[m.length];
        return recursiveGoodSpecial2(m,k,m2sum,m.length-1,best, maxDis);
    }
    public static int[] recursiveGoodSpecial2(int[] m, int k, int[] m2sum, int level, int[] best, int[] maxDis){
        int[] n = best;
        int min = Integer.MAX_VALUE;
        if(level ==0){
            best[0] = k;
        }else{
//            int start = Math.max(k-m2sum[level],0)/m[level];
            int start = Math.max(k-maxDis[level-1],0)/m[level];
            int end = k/m[level];
            if(level<m.length-1){
                end= Math.min(end,m[level+1]-1);
            }
            for (int i = end; i >= start; i--) {
                if(k-i*m[level]>maxDis[level-1]){
                    System.out.println("No");
                    return n;
                }
                best = recursiveGoodSpecial2(m,k-i*m[level], m2sum, level-1, best, maxDis);
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
    public static int[] recursiveGoodSpecial(int[] m, int k){
        int[] m2sum = new int[m.length];
        for (int i = 1; i < m2sum.length; i++) {
            m2sum[i]=m2sum[i-1]+m[i-1]*m[i];
        }
        int[] best = new int[m.length];
        return recursiveGoodSpecial(m,k,m2sum,m.length-1,best);
    }
    public static int[] recursiveGoodSpecial(int[] m, int k, int[] m2sum, int level, int[] best){
        int[] n = best;
        int min = Integer.MAX_VALUE;
        if(level ==0){
            best[0] = k;
        }else{
            int start = Math.max(k-m2sum[level],0)/m[level];
            int end = k/m[level];
            if(level<m.length-1){
                end= Math.min(end,m[level+1]);
            }
//            if(level==m.length-1)
//            System.out.println("Level "+level + " End " + end);
            for (int i = start; i < end+1; i++) {
//                if(i==251752) {
//                    System.out.println("TTT");
//                }
                best = recursiveGoodSpecial(m,k-i*m[level], m2sum, level-1, best);
                int sum=i;
                for (int j = 0; j < level; j++) {
                    sum += best[j];
                }
//                if(level==1&&i==2){
//                    System.out.println("ASD");
//                }
                if(Arrays.equals(new int[]{0, 2, 0, 0, 251752},best)) System.out.println("HHH");
                if(sum<min){
                    best[level]=i;
                    int summ=0;
                    for (int j = 0; j <= level; j++) {
                        summ+=m[j]*best[j];
                    }
                    if(summ!=k) continue;
                    min = sum;
                    n= best;
                    best = new int[m.length];
                }
            }
        }
        return n;
    }
    public static int[] recursiveGoodWrapped(int[] m, int k){
        int[] m2sum = new int[m.length];
        for (int i = 1; i < m2sum.length; i++) {
            m2sum[i]=m2sum[i-1]+m[i-1]*m[i];
        }
        int[] best = new int[m.length];
        return recursiveGoodWrapped(m,k,m2sum,m.length-1,best, new IntWrapper(Integer.MAX_VALUE));
    }
    public static int[] recursiveGood(int[] m, int k){
        int[] m2sum = new int[m.length];
        for (int i = 1; i < m2sum.length; i++) {
            m2sum[i]=m2sum[i-1]+m[i-1]*m[i];
        }
        int[] best = new int[m.length];
        return recursiveGood(m,k,m2sum,m.length-1,best);
    }
    public static int[] recursiveGood(int[] m, int k, int[] m2sum, int level, int[] best){
        int[] n = best;
        int min = Integer.MAX_VALUE;
        if(level ==0){
//            int fit = k/m[0];
//            if(fit*m[0]==k) {
//                best[0]=fit;
//            }
//            else best[0] = -1;
            best[0]=k;
        }else{
            int end = k/m[level];
            int start = Math.max(k-m2sum[level],0)/m[level];
            for (int i = start; i < end+1; i++) {
                best = recursiveGood(m,k-i*m[level], m2sum, level-1, best);
                if(best[0]==-1) continue;
                int sum=i;
                for (int j = 0; j < level; j++) {
                    sum += best[j];
                }
                if(sum<min){
                    min = sum;
                    n= best;
                    n[level]=i;
                    best = new int[m.length];
                }
            }
        }
        return n;
    }
    public static int[] recursiveGood2(int[] m, int k){
        int[] m2sum = new int[m.length];
        for (int i = 1; i < m2sum.length; i++) {
            m2sum[i]=m2sum[i-1]+m[i-1]*m[i];

        }
        int[] best = new int[m.length];
        return recursiveGood2(m,k,m2sum,m.length-1);
    }
    public static int[] recursiveGood2(int[] m, int k, int[] m2sum, int level){
        int[] n = new int[m.length];
        int min=Integer.MAX_VALUE;
        if(level ==0){
            int fit = k/m[0];
            if(fit*m[0]==k) {
                n[0]=fit;
            }
            else n[0] = -1;
        }else{
            int last = level;
            int loopCount = k/m[last];
            int start = Math.max(k-m2sum[level],0)/m[last];
            for (int i = start; i < loopCount+1; i++) {
                int[] newN =recursiveGood2(m,k-i*m[last], m2sum, level-1);
                if(newN[0]==-1) continue;
                int sum = i;
                for (int j = 0; j < level; j++) {
                    sum += newN[j];
                }
                if(sum<min){
                    min = sum;
//                    n= Arrays.copyOf(newN,newN.length+1);
                    n = newN;
                    n[last]=i;
                }
            }
        }
        return n;
    }
    private static class IntWrapper{
        int value;
        public IntWrapper(){}
        public IntWrapper(int value){this.value=value;}
    }
    public static int[] recursiveGoodWrapped(int[] m, int k, int[] m2sum, int level, int[] best, IntWrapper min){
        int[] n = best;
        if(level ==0){
            int fit = k/m[0];
            if(fit*m[0]==k) {
                best[0]=fit;
                min.value = fit;
            }
            else best[0] = -1;
        }else{
            int end = k/m[level];
            int start = Math.max(k-m2sum[level],0)/m[level];
            for (int i = start; i < end+1; i++) {
                IntWrapper sum=new IntWrapper(Integer.MAX_VALUE);
                best = recursiveGoodWrapped(m,k-i*m[level], m2sum, level-1, best, sum);
                if(best[0]==-1) continue;
                sum.value+=i;
                if(sum.value<min.value){
                    min.value = sum.value;
                    n= best;
                    n[level]=i;
                    best = new int[m.length];
                }
            }
        }
        return n;
    }
    public static int[] recursive(int[] m, int k){
        int[] n=new int[m.length];
        int min=Integer.MAX_VALUE;
        if(m.length==1){
            int fit = k/m[0];
            if(fit*m[0]==k) {
                n[0]=fit;
            }
            else n[0] = -1;
        }else{
            int last = m.length-1;
            int[] newM= Arrays.copyOf(m,last);
            int loopCount = k/m[last];
            int k2 = k;
            for (int i = 0; i < last; i++) {
                k2 -= m[i+1]*m[i];
            }
            if(k2<0)k2=0;
            int start = k2/m[last];

            for (int i = start; i < loopCount+1; i++) {
                if(k-i*m[last] <0) System.out.println(k-i*m[last]);
                int[] newN =recursive(newM,k-i*m[last]);
                if(newN[0]==-1) continue;
                int sum = Arrays.stream(newN).sum() + i;
                if(sum<min){
                    min = sum;
                    n = Arrays.copyOf(newN,newN.length+1);
                    n[last]=i;
                }
            }
        }
        return n;
    }
    public static int[] recursiveBruteforce(int[] m, int k){
        int[] n=new int[m.length];
        int min=Integer.MAX_VALUE;
        if(m.length==1){
            int fit = k/m[0];
            if(fit*m[0]==k) {
                n[0]=fit;
            }
            else n[0] = -1;
        }else{
            int last = m.length-1;
            int[] newM= Arrays.copyOf(m,last);
            int loopCount = k/m[last];
            for (int i = 0; i < loopCount+1; i++) {
                int[] newN =recursiveBruteforce(newM,k-i*m[last]);
                if(newN[0]==-1) continue;
                int sum = Arrays.stream(newN).sum() + i;
                if(sum<min){
                    min = sum;
                    n = Arrays.copyOf(newN,newN.length+1);
                    n[last]=i;
                }
            }
        }
        return n;
    }
    public static int[] bruteforceBad(int[]m, int k){
        int[]n=new int[m.length+1];
        int[]minArray=n;
        int min=Integer.MAX_VALUE;
        while(n[m.length]==0){
            int sum = 0;
            for (int i = 0; i < m.length; i++) {
                sum += n[i]*m[i];
            }
            if(sum==k){
                int nsum = 0;
                for (int i = 0; i < m.length; i++){
                    nsum+=n[i];
                }
                if(nsum<min){
                    min=nsum;
                    minArray= Arrays.copyOf(n,m.length);
                }
            }

            n[0]+=1;
            for (int i = 0; i < m.length; i++) {
                if(n[i]>=k/m[i]+1){
                    n[i+1]+=1;
                    n[i]=0;
                }
            }
        }
//        for (int i = 0; i < m.length; i++) {
//            for (int j = 0; j < k; j++) {
//                m[i]+=1;
//            }
//        }
        return minArray;
    }
}
